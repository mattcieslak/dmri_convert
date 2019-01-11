"""Handle merging and spliting of DSI files."""
import numpy as np
import os
import os.path as op
import nibabel as nb
from dipy.core.geometry import cart2sphere
from dipy.core.sphere import HemiSphere
from dipy.direction import peak_directions
import subprocess
from scipy.io.matlab import loadmat, savemat
from tqdm import tqdm
ODF_COLS = 20000  # Number of columns in DSI Studio odf split


def popen_run(arg_list):
    cmd = subprocess.Popen(arg_list, stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
    out, err = cmd.communicate()
    print(out)
    print(err)


def get_dsi_studio_ODF_geometry(odf_key):
    m = loadmat("odfs.mat")
    odf_vertices = m[odf_key + "_vertices"].T
    odf_faces = m[odf_key + "_faces"].T
    return odf_vertices, odf_faces


def mrtrix_to_dsistudio(mif_file, mask_file, output_fib, n_fibers=3):
    """Convert a MRTrix ``mif`` file containing sh coefficients to a DSI Studio fib file.

    This function uses ``sh2amp`` to get amplitude values for each direction in
    DSI Studio's ``odf8`` direction set. These values are masked and loaded into
    the "odfN" matrices in the fib file. The peak directions and peak indices are
    found using Dipy. **NOTE** the input data and mask file must be in LPS+
    orientation or none of this will work.

    Parameters:
    ============

    mif_file : str
        Path to a MRTrix file created by ``dwi2fod``. Contains sh coefficients
    mask_file : str
        Path to a nifti file that masks ``mif_file``
    output_fib : str
        Name of the output fib file. Must end with "fib", not "fib.gz"
    n_fibers : int
        The maximum number ODF maxima to extract per voxel

    """
    verts, faces = get_dsi_studio_ODF_geometry("odf8")
    num_dirs, _ = verts.shape
    hemisphere = num_dirs // 2
    x, y, z = verts[:hemisphere].T

    hs = HemiSphere(x=x, y=y, z=z)

    # Convert to DSI Studio LPS+ from MRTRIX3 RAS+
    _, theta, phi = cart2sphere(-x, -y, z)
    dirs_txt = op.join(os.getcwd(), "directions.txt")
    np.savetxt(dirs_txt, np.column_stack([phi, theta]))

    odf_amplitudes_nii = op.join(os.getcwd(), "amplitudes.nii.gz")
    popen_run(["sh2amp", "-force", "-nonnegative", mif_file, dirs_txt, odf_amplitudes_nii])

    if not op.exists(odf_amplitudes_nii):
        raise FileNotFoundError("Unable to create %s", odf_amplitudes_nii)
    amplitudes_img = nb.load(odf_amplitudes_nii)
    ampl_data = amplitudes_img.get_fdata()

    mask_img = nb.load(mask_file)

    if not np.allclose(mask_img.affine, amplitudes_img.affine):
        raise ValueError("Differing orientation between mask and amplitudes")
    if not mask_img.shape == amplitudes_img.shape[:3]:
        raise ValueError("Differing grid between mask and amplitudes")

    # Get the flat mask
    flat_mask = mask_img.get_data().flatten(order="F") > 0
    odf_array = ampl_data.reshape(-1, ampl_data.shape[3], order="F")
    masked_odfs = odf_array[flat_mask, :]
    n_odfs = masked_odfs.shape[0]
    peak_indices = np.zeros((n_odfs, n_fibers))
    peak_vals = np.zeros((n_odfs, n_fibers))

    dsi_mat = {}
    # Create matfile that can be read by dsi Studio
    dsi_mat['dimension'] = np.array(amplitudes_img.shape[:3])
    dsi_mat['voxel_size'] = np.array(amplitudes_img.header.get_zooms()[:3])
    n_voxels = int(np.prod(dsi_mat['dimension']))

    for odfnum in tqdm(range(masked_odfs.shape[0])):
        dirs, vals, indices = peak_directions(masked_odfs[odfnum], hs)
        for dirnum, (val, idx) in enumerate(zip(vals, indices)):
            if dirnum == n_fibers:
                break
            peak_indices[odfnum, dirnum] = idx
            peak_vals[odfnum, dirnum] = val

    for nfib in range(n_fibers):
        # fill in the "fa" values
        fa_n = np.zeros(n_voxels)
        fa_n[flat_mask] = peak_vals[:, nfib]
        dsi_mat['fa%d' % nfib] = fa_n.astype(np.float32)

        # Fill in the index values
        index_n = np.zeros(n_voxels)
        index_n[flat_mask] = peak_indices[:, nfib]
        dsi_mat['index%d' % nfib] = index_n.astype(np.int16)

    # Add in the ODFs
    num_odf_matrices = n_odfs // ODF_COLS
    split_indices = (np.arange(num_odf_matrices) + 1) * ODF_COLS
    odf_splits = np.array_split(masked_odfs, split_indices, axis=0)
    for splitnum, odfs in enumerate(odf_splits):
        dsi_mat['odf%d' % splitnum] = odfs.T.astype(np.float32)

    dsi_mat['odf_vertices'] = verts.T
    dsi_mat['odf_faces'] = faces.T
    dsi_mat['z0'] = np.array([1.])
    savemat("from_mrtrix.fib", dsi_mat, format='4', appendmat=False)
    os.remove("amplitudes.nii.gz")
    os.remove(dirs_txt)
