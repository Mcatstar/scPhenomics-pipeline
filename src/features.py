from pathlib import Path
import re

import numpy as np
from numpy.typing import NDArray
import pandas as pd
from scipy import ndimage as nd 
from tifffile import imread, imwrite
from cellpose import models
from loguru import logger
from tqdm import tqdm
import typer

from src.utils import object_quant, roi_quant
from src.config import PROCESSED_DATA_DIR, IMAGES_DATA_DIR

app = typer.Typer()


@app.command()
def main(
    # ---- REPLACE DEFAULT PATHS AS APPROPRIATE ----
    input_path: Path = IMAGES_DATA_DIR / "PLIN5 (1-424)", # 原始图像
    predictions_path: Path = PROCESSED_DATA_DIR / "PLIN5 (1-424)", # 模型推理输出数据
    output_path: Path = PROCESSED_DATA_DIR / "PLIN5 (1-424)",
    # -----------------------------------------
):
    # ---- REPLACE THIS WITH YOUR OWN CODE ----
    logger.info("Generating features from dataset...")
    logger.info("Pocessing predicted dataset...")
    # Organelle overlap
    tmp_list: list[Path] = list(Path(input_path).iterdir())
    img_list: list[Path] = []
    for tmp in tmp_list:
        if re.search('.tif', str(tmp)):
            img_list.append(tmp)

    logger.info("images found:\n"+"\n".join(map(lambda x: str(x.name), img_list)))
    logger.info(f"Number of images found: {len(img_list)}")

    for img_index in tqdm(range(len(img_list)), total=len(img_list), desc="CellPose Segmentation"):
        # Select image and load it ----------------------------------------------------
        img_path = img_list[img_index]
        img = imread(img_path)
        # Split channels ---------------------------------------------------------------
        actin_ch_number = int("3") - 1
        mito_ch_number = int("1") - 1
        lipid_ch_number = int("2") - 1
        actin_img = img[actin_ch_number, :, :]
        mito_img = img[mito_ch_number, :, :]
        lipid_img = img[lipid_ch_number, :, :]
        # load masks ---------------------------------------------------------------
        cell_mask_path = predictions_path / (img_list[img_index].name.replace('.tif', '_cell.tif'))
        cell_mask = imread(cell_mask_path)
        mito_mask_path = predictions_path / (img_list[img_index].name.replace('.tif', '_mito.tif'))
        mito_mask = imread(mito_mask_path)
        lipid_mask_path = predictions_path / (img_list[img_index].name.replace('.tif', '_lipid.tif'))
        lipid_mask = imread(lipid_mask_path)
        cv_mask_path = predictions_path / (img_list[img_index].name.replace('.tif', '_cv_mask.tif'))
        cv_mask = imread(cv_mask_path)
        pv_mask_path = predictions_path / (img_list[img_index].name.replace('.tif', '_pv_mask.tif'))
        pv_mask = imread(pv_mask_path)
        mito_edt_path = predictions_path / (img_list[img_index].name.replace('.tif', '_mito_edt.tif'))
        mito_edt = imread(mito_edt_path)
        lipid_edt_path = predictions_path / (img_list[img_index].name.replace('.tif', '_lipid_edt.tif'))
        lipid_edt = imread(lipid_edt_path)

        # Organelle overlap----------------------------------------
        mito_binary = np.where(mito_mask > 0, 1, 0)
        lipid_binary = np.where(lipid_mask > 0 , 1, 0)
        overlap_tmp = np.add(mito_binary, lipid_binary)
        overlap_mask = np.where(overlap_tmp == 2, cell_mask, 0)

        overlap_mask_path = predictions_path / (img_list[img_index].name.replace('.tif', '_overlap.tif'))
        imwrite(overlap_mask_path, overlap_mask)

        # EDT maps--------------------------------------------------
        cv_mask = np.where(cv_mask == 1, 0, 1)
        cv_edt = nd.distance_transform_edt(cv_mask)

        pv_mask = np.where(pv_mask == 1, 0, 1)
        pv_edt = nd.distance_transform_edt(pv_mask)

        mito_binary = np.where(mito_mask == 0, 1, 0)
        mito_edt = nd.distance_transform_edt(mito_binary)

        lipid_binary = np.where(lipid_mask == 0, 1, 0)
        lipid_edt = nd.distance_transform_edt(lipid_binary)

        # Quantification----------------------------------------------------
        pv_df_tmp = roi_quant(pv_mask, cell_mask, mito_img, lipid_img, cv_edt, pv_edt, mito_edt)
        pv_cv_dist = pv_df_tmp['cv_max_dist'][0]
        logger.info(f"pv_cv_dist: {pv_cv_dist}")

        cell_df_tmp = roi_quant(cell_mask, cell_mask, mito_img, lipid_img, cv_edt, pv_edt, mito_edt)
        cell_df = pd.DataFrame()
        cell_df['cell_id'] = cell_df_tmp['roi_id']
        cell_df['cell_area'] = cell_df_tmp['area']
        cell_df['cell_centroid_x'] = cell_df_tmp['centroid_x']
        cell_df['cell_centroid_y'] = cell_df_tmp['centroid_y']
        cell_df['cell_mito_int'] = cell_df_tmp['mito_int']
        cell_df['cell_lipid_int'] = cell_df_tmp['lipid_int']
        cell_df['cell_cv_min_dist'] = cell_df_tmp['cv_min_dist']
        cell_df['cell_cv_max_dist'] = cell_df_tmp['cv_max_dist']
        cell_df['cell_cv_mean_dist'] = cell_df_tmp['cv_mean_dist']
        cell_df['cell_pv_min_dist'] = cell_df_tmp['pv_min_dist']
        cell_df['cell_pv_max_dist'] = cell_df_tmp['pv_max_dist']
        cell_df['cell_pv_mean_dist'] = cell_df_tmp['pv_mean_dist']

        logger.info(f"cell_df.head(): {cell_df.head()}")
        logger.info(f"cell_df.shape: {cell_df.shape}")

        mito_df_tmp = roi_quant(mito_mask, cell_mask, mito_img, lipid_img, cv_edt, pv_edt, lipid_edt)
        mito_df = mito_df_tmp.merge(cell_df, how = 'outer', on = 'cell_id')
        mito_df['organelle'] = 'mito'
        logger.info(f"mito_df.head(): {mito_df.head()}")
        logger.info(f"mito_df.shape: {mito_df.shape}")

        lipid_df_tmp = roi_quant(lipid_mask, cell_mask, mito_img, lipid_img, cv_edt, pv_edt, mito_edt)
        lipid_df = lipid_df_tmp.merge(cell_df, how = 'outer', on = 'cell_id')
        lipid_df['organelle'] = 'lipid'
        logger.info(f"lipid_df.head(): {lipid_df.head()}")
        logger.info(f"lipid_df.shape: {lipid_df.shape}")

        overlap_df_tmp = roi_quant(overlap_mask, cell_mask, mito_img, lipid_img, cv_edt, pv_edt, mito_edt)
        overlap_df = pd.DataFrame() 
        overlap_df['cell_id'] = overlap_df_tmp['roi_id']
        overlap_df['overlap_area'] = overlap_df_tmp['area']
        logger.info(f"overlap_df.head(): {overlap_df.head()}")
        logger.info(f"overlap_df.shape: {overlap_df.shape}")

        df = pd.concat([mito_df, lipid_df])
        df = df.merge(overlap_df, how = 'outer', on = 'cell_id')
        df['pv_cv_dist'] = pv_cv_dist
        logger.info(f"head(5) of Combine dataframes: {df.head()}")
        logger.info(f"shape of Combine dataframes: {df.shape}")

        df_path = output_path / (img_list[img_index].name.replace('.tif', '_output.csv'))
        df.to_csv(df_path, header = True)
        logger.success("Statistics csv generation complete.")
    # -----------------------------------------


if __name__ == "__main__":
    app()
