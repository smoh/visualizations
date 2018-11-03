# %matplotlib inline
# %config InlineBackend.figure_format = 'retina'

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
import palettable

from astropy.table import Table
import astropy.units as u
import astropy.coordinates as coord
from astropy.io import fits

# %% Data: combine APOGEE with spectrophotometric parallaxes
t = fits.open("dr15/allStar-t9-l31c-58158.fits")[1].data
specd = fits.getdata('data_HoggEilersRix2018.fits')
print('apogee has {:d} rows'.format(len(t)))
print('specd has {:d} rows'.format(len(specd)))

merged = pd.merge(
    pd.DataFrame({'TMASS_ID':t['APOGEE_ID']}).reset_index(),
    pd.DataFrame({'TMASS_ID':specd['2MASS_ID']}).reset_index(),
    on='TMASS_ID')
print('N match =', len(merged))

def fits_rec_to_dataframe(fits_rec):
    return pd.DataFrame(fits_rec.tolist(), columns=fits_rec.dtype.names)

df = pd.concat(
    [fits_rec_to_dataframe(t[merged.index_x]),
     fits_rec_to_dataframe(specd[merged.index_y])],
    axis=1)

# %%
coords = coord.ICRS(
    df['RA'].values*u.deg, df['DEC'].values*u.deg,
    1*u.kpc/df['spec_parallax'].values)
xyz = coords.transform_to(coord.Galactic).cartesian.xyz

logg_cut = (1 < df['LOGG']) & (df['LOGG'] < 3.8)
teff_cut = (3500 < df['TEFF']) & (df['TEFF'] < 5500)
snr_cut = df['SNR'] > 80
bad_flag = ~(df['ASPCAPFLAG'] & 2**23).astype(bool)
target1_flag = (
    (df['APOGEE_TARGET1'] & 2**11)
    | (df['APOGEE_TARGET1'] & 2**12)
    | (df['APOGEE_TARGET1'] & 2**13)).astype(bool)
has_feh = (df['FE_H']!=-9999) | (df['FE_H_ERR']!=-9999)
print(logg_cut.sum(), teff_cut.sum(), snr_cut.sum())
print((logg_cut & teff_cut & snr_cut).sum())
print(target1_flag.sum())
print((snr_cut & bad_flag).sum())


# %%
# fig, ax = plt.subplots(1,1)
# ax.hist2d(t['TEFF'][merged.index_x], t['LOGG'][merged.index_x],
#     bins=(np.linspace(3000, 6000, 64), np.linspace(-1, 3.8, 64)),
#     norm=colors.LogNorm());
#
#
#
# fig, ax = plt.subplots(1, 2, figsize=(10, 4))
# ax[0].hist2d(xyz[0], xyz[1], bins=128, norm=colors.LogNorm())
# ax[1].hist2d(xyz[0], xyz[2], bins=128, norm=colors.LogNorm());
#
# print((snr_cut[merged.index_x] & bad_flag[merged.index_x]).sum())
#
# import seaborn as sns
# from scipy import stats
#
#
# s, xedges, yedges, binnumber = stats.binned_statistic_2d(
#     xyz[0], xyz[1], tt['FE_H'], bins=128, statistic='median')
#
#
# plt.pcolormesh(xedges[:-1], yedges[:-1], s.T, vmin=-2);
#
# plt.rc('figure', dpi=120)
# plt.figure()
# plt.imshow(s.T, vmin=-2, origin='lower',)
# plt.colorbar();
#
#
# s, xedges, yedges, binnumber = stats.binned_statistic_2d(
#     xyz[0], xyz[1], tt['AL_H'], bins=128, statistic='median')
#
# plt.rc('figure', dpi=120)
# plt.figure()
# plt.imshow(s.T, vmin=-2)
# plt.colorbar();
#
# plt.scatter(tt['FE_H'], tt['FE_H_ERR'], s=4, alpha=.1);
# plt.ylim(0,.1)
# plt.xlim(-3,1)
