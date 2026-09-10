# Import libraries---------------------------------------------------------
from typing import Any

import numpy as np 
from numpy.typing import NDArray
import pandas as pd 
from skimage.measure import regionprops
from tqdm import tqdm
from cellpose import models
from loguru import logger


# Define functions---------------------------------------------------------
def model_apply(model: models.CellposeModel, img: NDArray) -> NDArray | Any:
    """Apply the Cellpose model to an image.
    You will get a tuple containing (masks, flows, styles) as output.
    The mask is the most important output,
    which contains the segmented objects in the image.
    """
    logger.info("Applying Cellpose model to image...")
    mask, _, _ = model.eval(img, channels = [0, 0])
    
    return mask


def object_quant(lbl_img, lipid_img) -> pd.DataFrame:
    """
    Quantify mean lipid intensity for each object in a labeled image.

    This function computes the mean intensity of the lipid channel for every
    segmented region (object) in the input label image. It is typically used
    as a precursor to filtering objects based on lipid signal intensity.

    Parameters
    ----------
    lbl_img : np.ndarray
        2D labeled image where each unique positive integer represents a
        distinct segmented object (e.g., lipid droplets).
    lipid_img : np.ndarray
        2D grayscale image corresponding to the lipid fluorescence channel.
        Must have the same spatial dimensions as lbl_img.

    Returns
    -------
    pd.DataFrame
        A DataFrame with the following columns:
        - 'roi_id' : int, the unique label of each object.
        - 'lipid_int' : float, the mean pixel intensity of the lipid channel
          within the object's region.
    """
    df = pd.DataFrame() 
    
    roi_id = [] 
    lipid_int = [] 
    
    roiprops = regionprops(lbl_img, intensity_image = lipid_img)
    
    for roi in tqdm(range(len(roiprops))):
        roi_id.append(roiprops[roi].label)
        lipid_int.append(roiprops[roi].intensity_mean)
        
    df['roi_id'] = roi_id
    df['lipid_int'] = lipid_int
    
    return df


def object_filter(
        lbl_img, int_df: pd.DataFrame, int_thresh: int
        ) -> NDArray:
    """
    Filter labeled objects based on a lipid intensity threshold.

    This function removes objects whose mean lipid intensity falls below a
    specified threshold by setting their pixel values to zero. It preserves
    the original labels for objects that pass the filter.

    Parameters
    ----------
    lbl_img : np.ndarray
        2D labeled image with integer labels for each object.
    int_df : pd.DataFrame
        DataFrame containing object intensity data. Must include a column
        named 'roi_id' (object labels) and a column named 'lipid_int'
        (mean lipid intensity for each object).
    int_thresh : int
        Intensity threshold. Objects with mean lipid intensity less than or
        equal to this value will be removed (set to 0).

    Returns
    -------
    np.ndarray
        A new labeled image of the same shape as lbl_img, where objects
        with lipid_int <= int_thresh have been set to 0. Objects that pass
        the threshold retain their original labels.
    """
    roi_id = lbl_img.reshape(-1)
    
    mask_df = pd.DataFrame()
    mask_df['roi_id'] = roi_id
    mask_df = pd.merge(mask_df, int_df, how = 'left', on = ['roi_id'])
    
    mask_tmp = mask_df['lipid_int'].to_numpy()
    mask_tmp = mask_tmp.reshape(lbl_img.shape)
    
    mask = np.where(mask_tmp <= int_thresh, 0, lbl_img)
    
    return mask


def roi_quant(
        lbl_img: NDArray, cell_img: NDArray, mito_img: NDArray, lipid_img: NDArray, 
        cv_edt, pv_edt, organelle_edt
        ) -> pd.DataFrame:
    """
    Extract comprehensive morphometric, intensity, and spatial features for
    each object in a labeled image.

    This function computes a wide range of quantitative features for each
    segmented object (e.g., a cell) using the input label image and
    corresponding fluorescence and distance transform maps.

    Parameters
    ----------
    lbl_img : np.ndarray
        2D labeled image where each unique positive integer represents a
        distinct segmented object (e.g., a cell).
    cell_img : np.ndarray
        2D image used for identifying cell boundaries (e.g., actin channel).
        Used to compute intensity features for each object.
    mito_img : np.ndarray
        2D fluorescence image of mitochondria. Used to compute mean
        mitochondrial intensity within each object.
    lipid_img : np.ndarray
        2D fluorescence image of lipid droplets. Used to compute mean
        lipid intensity within each object.
    cv_edt : np.ndarray
        2D Euclidean Distance Transform (EDT) map from the central vein (CV).
        Pixel values represent distance to the nearest CV boundary.
    pv_edt : np.ndarray
        2D Euclidean Distance Transform (EDT) map from the portal vein (PV).
        Pixel values represent distance to the nearest PV boundary.
    organelle_edt : np.ndarray
        2D Euclidean Distance Transform (EDT) map from an organelle of
        interest (e.g., mitochondria or lipid droplets). Pixel values
        represent distance to the nearest organelle surface.

    Returns
    -------
    pd.DataFrame
        A DataFrame with one row per object (ROI) containing the following
        feature groups:

        **Identification and Morphometrics:**
        - 'roi_id' : int, unique label of the object.
        - 'cell_id' : int, rounded mean intensity from cell_img (used as
          a secondary identifier).
        - 'area' : int, number of pixels in the object.
        - 'perimeter' : float, perimeter length of the object.
        - 'eccentricity' : float, eccentricity of the fitted ellipse
          (0 = circle, closer to 1 = elongated).
        - 'solidity' : float, ratio of area to convex hull area
          (1 = perfectly convex).
        - 'axis_major_length' : float, major axis length of fitted ellipse.
        - 'axis_minor_length' : float, minor axis length of fitted ellipse.
        - 'feret_diameter_max' : float, maximum Feret diameter (caliper
          distance).
        - 'centroid_x' : float, X-coordinate (column) of the object centroid.
        - 'centroid_y' : float, Y-coordinate (row) of the object centroid.

        **Intensity Features:**
        - 'mito_int' : float, mean mitochondrial fluorescence intensity.
        - 'lipid_int' : float, mean lipid fluorescence intensity.

        **Spatial Features (Central Vein):**
        - 'cv_min_dist' : float, minimum distance from object to CV.
        - 'cv_max_dist' : float, maximum distance from object to CV.
        - 'cv_mean_dist' : float, mean distance from object to CV.

        **Spatial Features (Portal Vein):**
        - 'pv_min_dist' : float, minimum distance from object to PV.
        - 'pv_max_dist' : float, maximum distance from object to PV.
        - 'pv_mean_dist' : float, mean distance from object to PV.

        **Spatial Features (Organelle Proximity):**
        - 'organelle_min_dist' : float, minimum distance from object to
          the nearest organelle.
        - 'organelle_max_dist' : float, maximum distance from object to
          the nearest organelle.
        - 'organelle_mean_dist' : float, mean distance from object to the
          nearest organelle.

    Notes
    -----
    This function uses scikit-image's `regionprops` internally and iterates
    over all objects. For large images with many objects, this operation
    can be time-consuming. Progress is reported via a tqdm progress bar.
    """
    df = pd.DataFrame() 
    
    roi_id = [] 
    cell_id = [] 
    area = [] 
    perimeter = []
    eccentricity = [] 
    solidity = [] 
    axis_major_length = []
    axis_minor_length = [] 
    feret_diameter_max = [] 
    centroid_x = [] 
    centroid_y = [] 
    mito_int = [] 
    lipid_int = [] 
    cv_min_dist = [] 
    cv_max_dist = [] 
    cv_mean_dist = [] 
    pv_min_dist = [] 
    pv_max_dist = [] 
    pv_mean_dist = []
    organelle_min_dist = []
    organelle_max_dist = [] 
    organelle_mean_dist = [] 
    
    cellprops = regionprops(lbl_img, intensity_image = cell_img)
    mitoprops = regionprops(lbl_img, intensity_image = mito_img)
    lipidprops = regionprops(lbl_img, intensity_image = lipid_img)
    cvprops = regionprops(lbl_img, intensity_image = cv_edt)
    pvprops = regionprops(lbl_img, intensity_image = pv_edt)
    organelleprops = regionprops(lbl_img, intensity_image = organelle_edt)
    
    for roi in tqdm(range(len(cellprops))):
        roi_id.append(cellprops[roi].label)
        cell_id.append(round(cellprops[roi].intensity_mean))
        area.append(cellprops[roi].area)
        perimeter.append(cellprops[roi].perimeter)
        eccentricity.append(cellprops[roi].eccentricity)
        solidity.append(cellprops[roi].solidity)
        axis_major_length.append(cellprops[roi].axis_major_length)
        axis_minor_length.append(cellprops[roi].axis_minor_length)
        feret_diameter_max.append(cellprops[roi].feret_diameter_max)
        centroid_x.append(cellprops[roi].centroid[1])
        centroid_y.append(cellprops[roi].centroid[0])
        mito_int.append(mitoprops[roi].intensity_mean)
        lipid_int.append(lipidprops[roi].intensity_mean)
        cv_min_dist.append(cvprops[roi].intensity_min)
        cv_max_dist.append(cvprops[roi].intensity_max)
        cv_mean_dist.append(cvprops[roi].intensity_mean)
        pv_min_dist.append(pvprops[roi].intensity_min)
        pv_max_dist.append(pvprops[roi].intensity_max)
        pv_mean_dist.append(pvprops[roi].intensity_mean)
        organelle_min_dist.append(organelleprops[roi].intensity_min)
        organelle_max_dist.append(organelleprops[roi].intensity_max)
        organelle_mean_dist.append(organelleprops[roi].intensity_mean)
        
    df['roi_id'] = roi_id
    df['cell_id'] = cell_id
    df['area'] = area
    df['perimeter'] = perimeter
    df['eccentricity'] = eccentricity
    df['solidity'] = solidity
    df['axis_major_length'] = axis_major_length
    df['axis_minor_length'] = axis_minor_length
    df['feret_diameter_max'] = feret_diameter_max
    df['centroid_x'] = centroid_x
    df['centroid_y'] = centroid_y
    df['mito_int'] = mito_int
    df['lipid_int'] = lipid_int
    df['cv_min_dist'] = cv_min_dist
    df['cv_max_dist'] = cv_max_dist
    df['cv_mean_dist'] = cv_mean_dist
    df['pv_min_dist'] = pv_min_dist
    df['pv_max_dist'] = pv_max_dist
    df['pv_mean_dist'] = pv_mean_dist
    df['organelle_min_dist'] = organelle_min_dist
    df['organelle_max_dist'] = organelle_max_dist
    df['organelle_mean_dist'] = organelle_mean_dist
    
    return df 