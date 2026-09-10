from pathlib import Path
import re

import numpy as np
from numpy.typing import NDArray
from scipy import ndimage as nd 
from tifffile import imread, imwrite
from cellpose import models
from loguru import logger
from tqdm import tqdm
import typer

from src.utils import model_apply, object_quant, object_filter
from src.config import MODELS_DIR, IMAGES_DATA_DIR, PROCESSED_DATA_DIR

app = typer.Typer()


class Segment:
    def __init__(self, model_path: Path, img: NDArray, gpu: bool = False):
        self.model_path = model_path
        self.img = img
        # Load your model here (e.g., using joblib or pickle)
        logger.info(f"Loading model from {self.model_path}...")
        self.model = models.CellposeModel(gpu = gpu, pretrained_model = str(self.model_path)) # pyright: ignore[reportArgumentType]

    def predict(self):
        # Perform prediction using the loaded model
        logger.info("Performing prediction...")
        mask = model_apply(self.model, self.img)  # get mask
        predictions = np.uint16(mask)
        return predictions # return mask

    def save_mask(self, mask: NDArray, output_path: Path):
        # Save the mask to the specified output path
        logger.info(f"Saving mask to {output_path}...")
        imwrite(output_path, mask)
        logger.success(f"Mask saved to {output_path}.")


class CellSegmenter(Segment):
    def __init__(self, model_path: Path, mito_img: NDArray, actin_img: NDArray, gpu: bool):
        super().__init__(model_path, mito_img, gpu)
        self.cell_model_path = model_path
        self.mito_img = mito_img
        self.actin_img = actin_img
        self.cell_model = models.CellposeModel(gpu = False, pretrained_model = str(self.cell_model_path))# pyright: ignore[reportArgumentType]
        # Additional initialization for CellSegmenter if needed

    def predict(self):
        logger.info("Performing cell segmentation prediction...")
        mito_tmp = nd.gaussian_filter(self.mito_img, sigma = 50)
        mito_tmp = mito_tmp * 3
        tmp_img = np.add(self.actin_img, mito_tmp)

        cell_mask = model_apply(self.cell_model, tmp_img)
        cell_mask = np.uint16(cell_mask)
        return cell_mask

    def save_mask(self, cell_mask, cell_mask_path: Path):
        # Save the mask to the specified output path
        logger.info(f"Saving mask to {cell_mask_path}...")
        imwrite(cell_mask_path, cell_mask)
        logger.success(f"Mask saved to {cell_mask_path}.")


class MitoSegmenter(Segment):
    def __init__(self, model_path: Path, mito_img: NDArray, gpu: bool = False):
        super().__init__(model_path, mito_img, gpu)
        self.mito_model_path = model_path
        self.mito_img = mito_img
        self.mito_model = models.CellposeModel(gpu = False, pretrained_model = str(self.mito_model_path))# pyright: ignore[reportArgumentType]
        # Additional initialization for CellSegmenter if needed

    def predict(self):
        logger.info("Performing mito segmentation prediction...")
        mito_mask = model_apply(self.mito_model, self.mito_img)
        mito_mask = np.uint16(mito_mask)
        return mito_mask

    def save_mask(self, mito_mask, mito_mask_path: Path):
        # Save the mask to the specified output path
        logger.info(f"Saving mask to {mito_mask_path}...")
        imwrite(mito_mask_path, mito_mask)
        logger.success(f"Mask saved to {mito_mask_path}.")


class LipidSegmenter(Segment):
    def __init__(self, model_path: Path, lipid_img: NDArray, gpu: bool = False):
        super().__init__(model_path, lipid_img, gpu)
        self.lipid_model_path = model_path
        self.lipid_img = lipid_img
        self.lipid_model = models.CellposeModel(gpu = False, pretrained_model = str(self.lipid_model_path))# pyright: ignore[reportArgumentType]
        # Additional initialization for CellSegmenter if needed

    def predict(self):
        logger.info("Performing lipid segmentation prediction...")
        lipid_mask_tmp = model_apply(self.lipid_model, self.lipid_img)
        lipid_mask_tmp = np.uint16(lipid_mask_tmp)

        lipid_df_tmp = object_quant(lipid_mask_tmp, self.lipid_img)
        lipid_mask = object_filter(lipid_mask_tmp, lipid_df_tmp, 10)
        return lipid_mask

    def save_mask(self, lipid_mask, lipid_mask_path: Path):
        # Save the mask to the specified output path
        logger.info(f"Saving mask to {lipid_mask_path}...")
        imwrite(lipid_mask_path, lipid_mask)
        logger.success(f"Mask saved to {lipid_mask_path}.")


class quant:
    def __init__(self, model_path: Path):
        self.model_path = model_path
        # Load your model here (e.g., using joblib or pickle)
        logger.info(f"Loading model from {self.model_path}...")
        # self.model = load_model(self.model_path)  # Placeholder for actual model loading

    def predict(self, features):
        # Perform prediction using the loaded model
        logger.info("Performing prediction...")
        # predictions = self.model.predict(features)  # Placeholder for actual prediction
        predictions = [0] * len(features)  # Dummy predictions for illustration
        return predictions


@app.command()
def main(
    # ---- REPLACE DEFAULT PATHS AS APPROPRIATE ----
    input_path: Path = IMAGES_DATA_DIR / "PLIN5 variants on CNTR diet" / "PLIN5 (1-424)",
    model_path: Path = MODELS_DIR,
    predictions_path: Path = PROCESSED_DATA_DIR / "PLIN5 (1-424)",
    # -----------------------------------------
):
    # ---- REPLACE THIS WITH YOUR OWN CODE ----
    logger.info("Performing inference for model...")
    logger.info(f"""Input path: {input_path}\nOutput path: {predictions_path}
    """)

    # Load image-----------------------------------------------------------------------------------------------
    tmp_list: list[Path] = list(Path(input_path).iterdir())

    img_list: list[Path] = []
    for tmp in tmp_list:
        if re.search('.tif', str(tmp)):
            img_list.append(tmp)
            
    logger.info("images found:\n"+"\n".join(map(lambda x: str(x.name), img_list)))
    logger.info(f"Number of images found: {len(img_list)}")
    # img_index = int(input('Select image: '))
    for img_index in tqdm(range(len(img_list)), total=len(img_list), desc="CellPose Segmentation"):
        # Select image and load it
        logger.info(f"Selected image: {img_list[img_index].name}")
        img_path: Path = img_list[img_index] # img_list have full path
        logger.info(f"Image path: {img_path}")
        img = imread(img_path)
        logger.info(img.shape)

        # Split channels ---------------------------------------------------------------
        # 正确分辨通道是模型正确运行的基础，FIJI分离通道观察形状以辅助区分：Image > Color > Split Channels
        actin_ch_number = int("3") - 1
        mito_ch_number = int("1") - 1
        lipid_ch_number = int("2") - 1
        actin_img = img[actin_ch_number, :, :]
        mito_img = img[mito_ch_number, :, :]
        lipid_img = img[lipid_ch_number, :, :]

        # Cell segmentation ---------------------------------------------------------------
        cell_segmenter = CellSegmenter(model_path / "brownl_cell_01", mito_img, actin_img, gpu = False)
        cell_mask = cell_segmenter.predict()
        cell_mask_path = predictions_path / (img_list[img_index].name.replace('.tif', '_cell.tif'))
        cell_segmenter.save_mask(cell_mask, cell_mask_path)

        # Mito segmentation ---------------------------------------------------------------
        mito_segmenter = MitoSegmenter(model_path / "brownl_mito_03", mito_img, gpu = False)
        mito_mask = mito_segmenter.predict()
        mito_mask_path = predictions_path / (img_list[img_index].name.replace('.tif', '_mito.tif'))
        mito_segmenter.save_mask(mito_mask, mito_mask_path)

        # Lipid segmentation ---------------------------------------------------------------
        lipid_segmenter = LipidSegmenter(model_path / "brownl_lipid_01", lipid_img, gpu = False)
        lipid_mask = lipid_segmenter.predict()
        lipid_mask_path = predictions_path / (img_list[img_index].name.replace('.tif', '_lipid.tif'))
        lipid_segmenter.save_mask(lipid_mask, lipid_mask_path)

    logger.success("Inference complete.")
    # -----------------------------------------


if __name__ == "__main__":
    app()
