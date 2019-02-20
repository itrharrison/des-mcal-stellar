from astropy.table import Table
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
import os

rc('text', usetex=True)
rc('font', family='serif')
rc('font', size=11)

plt.close('all') # tidy up any unshown plots

def make_shape_cuts(cat, snr=[10,100], flags=0, size_ratio=0.5):
  flag_cut = cat['flags']==flags
  snr_cut = np.less(cat['snr'], snr[1])*np.greater(cat['snr'], snr[0])
  size_cut = np.greater(cat['T']/cat['mcal_psf_T'], size_ratio)
  shape_cut = flag_cut*snr_cut*size_cut

  return shape_cut

hsc_fname = 'hsc-sxds.fits'
gold_fname = 'sxds.y3_gold_2_2.fits'
mcal_fname = 'sxds_metacal_unsheared_Y3_mastercat_v2_6_20_18_subsampled.fits'

hsc_gold_match_fname = 'sxds.hsc_gold_matches.fits'
mcal_gold_match_fname = 'sxds.mcal_gold_matches.fits'
mcal_hsc_gold_match_fname = 'sxds.hsc_mcal_gold_matches.fits'

hsc = Table.read(hsc_fname)
gold = Table.read(gold_fname)
mcal = Table.read(mcal_fname)

stilts_cmd = 'java -jar /Applications/TOPCAT.app/Contents/Resources/Java/topcat-full.jar -stilts'
stilts_sky_join_cmd = ' tmatch2 in1={0} in2={1} out={2} matcher=sky values1="ra dec" values2="ra dec" params="0.5"'
stilts_id_join_cmd = ' tmatch2 in1={0} in2={1} out={2} matcher=exact values1="coadd_object_id" values2="COADD_OBJECT_ID"'


# cross match all the catalogues
if not os.path.exists(hsc_gold_match_fname):
  os.system(stilts_cmd+stilts_sky_join_cmd.format(hsc_fname, gold_fname, hsc_gold_match_fname))

hsc_gold_matches = Table.read(hsc_gold_match_fname)

if not os.path.exists(mcal_gold_match_fname):
  os.system(stilts_cmd+stilts_id_join_cmd.format(mcal_fname, gold_fname, mcal_gold_match_fname))

mcal_gold_matches = Table.read(mcal_gold_match_fname)

if not os.path.exists(mcal_hsc_gold_match_fname):
  os.system(stilts_cmd+stilts_id_join_cmd.format(mcal_fname, hsc_gold_match_fname, mcal_hsc_gold_match_fname))

mcal_hsc_gold_matches = Table.read(mcal_hsc_gold_match_fname)

# shape cuts to catalogues
shape_mcal_gold_cut = make_shape_cuts(mcal_gold_matches)
shape_mcal_hsc_gold_cut = make_shape_cuts(mcal_hsc_gold_matches)

# y1 spreadmodel cuts to mcal-gold
star_spread_cut = np.less(np.abs(mcal_gold_matches['SPREAD_MODEL_I'] + (5./3.)*mcal_gold_matches['SPREADERR_MODEL_I']), 0.002)
gal_spread_cut = ~star_spread_cut

# hsc star cuts to mcal-hsc-gold
# need to recalculate shape cuts!
star_hsc_cut = mcal_hsc_gold_matches['iclassification_extendedness']==0
gal_hsc_cut = ~star_hsc_cut
R11_hsc_cut = np.less(mcal_hsc_gold_matches['R11'], 3.)*np.greater(mcal_hsc_gold_matches['R11'], -3)

# need to think about completeness across this
N_mcal = len(mcal)
N_gold = len(gold)
N_hsc = len(hsc)

N_mcal_gold = len(mcal_gold_matches)
N_hsc_gold = len(hsc_gold_matches)
N_mcal_hsc_gold = len(mcal_hsc_gold_matches)

N_shape = np.sum(shape_mcal_gold_cut)
N_shape_hsc = np.sum(shape_mcal_hsc_gold_cut)
N_hsc_shape_stars = np.sum(star_hsc_cut*shape_mcal_hsc_gold_cut)
N_spread_shape_stars = np.sum(star_spread_cut*shape_mcal_gold_cut)

print(N_mcal, N_gold, N_hsc)
print(N_mcal_gold, N_hsc_gold, N_mcal_hsc_gold)

print(N_shape, N_hsc_shape_stars, N_spread_shape_stars)
print(N_shape, N_hsc_shape_stars/N_shape, N_spread_shape_stars/N_shape)

print('HSC contamination: {0:.3f}\\%'.format(100*N_hsc_shape_stars/N_shape))
print('Y1Spreadmodel contamination: {0:.3f}\\%'.format(100*N_spread_shape_stars/N_shape))

plt.figure(1, figsize=(2*4.5, 3.75))
plt.subplot(121)
plt.title('Galaxies')
plt.hist2d(mcal_gold_matches['MAG_AUTO_I'][gal_spread_cut], mcal_gold_matches['SPREAD_MODEL_I'][gal_spread_cut], range=[[18, 25],[-0.02, 0.05]], bins=50, cmap='Blues')
plt.xlabel('MAG AUTO I')
plt.ylabel('SPREAD MODEL I')
plt.subplot(122)
plt.title('Stars')
plt.hist2d(mcal_gold_matches['MAG_AUTO_I'][star_spread_cut], mcal_gold_matches['SPREAD_MODEL_I'][star_spread_cut], range=[[18, 25],[-0.02, 0.05]], bins=50, cmap='Oranges')
plt.xlabel('MAG AUTO I')
#plt.ylabel('SPREAD MODEL I')
plt.suptitle('Y1 Spread Model')
plt.savefig('plots/stargal-y1spreadmodel.png', dpi=300, bbox_inches='tight')

plt.figure(2, figsize=(2*4.5, 3.75))
plt.subplot(121)
plt.title('Galaxies')
plt.hist2d(mcal_hsc_gold_matches['MAG_AUTO_I'][gal_hsc_cut], mcal_hsc_gold_matches['SPREAD_MODEL_I'][gal_hsc_cut], range=[[18, 25],[-0.02, 0.05]], bins=50, cmap='Blues')
plt.xlabel('MAG AUTO I')
plt.ylabel('SPREAD MODEL I')
plt.subplot(122)
plt.title('Stars')
plt.hist2d(mcal_hsc_gold_matches['MAG_AUTO_I'][star_hsc_cut], mcal_hsc_gold_matches['SPREAD_MODEL_I'][star_hsc_cut], range=[[18, 25],[-0.02, 0.05]], bins=50, cmap='Oranges')
plt.xlabel('MAG AUTO I')
#plt.ylabel('SPREAD MODEL I')
plt.suptitle('HSC Extendedness')
plt.savefig('plots/stargal-hsc.png', dpi=300, bbox_inches='tight')

R11_gold_cut = np.less(mcal_gold_matches['R11'], 3.)*np.greater(mcal_gold_matches['R11'], -3)*(mcal_gold_matches['flags']==0)
R11_hsc_cut = np.less(mcal_hsc_gold_matches['R11'], 3.)*np.greater(mcal_hsc_gold_matches['R11'], -3)*(mcal_hsc_gold_matches['flags']==0)

plt.figure(7, figsize=(2*4.5, 3.75))
plt.subplot(121)
plt.hist(mcal_gold_matches['R11'][R11_gold_cut*shape_mcal_gold_cut], bins=50, histtype='step', normed=True, label='MCAL Shape Objects')
plt.hist(mcal_hsc_gold_matches['R11'][R11_hsc_cut*shape_mcal_hsc_gold_cut*star_hsc_cut], bins=50, histtype='step', normed=True, label='MCAL Shape, HSC Stars')
plt.hist(mcal_hsc_gold_matches['R11'][R11_hsc_cut*shape_mcal_hsc_gold_cut*gal_hsc_cut], bins=50, histtype='step', normed=True, label='MCAL Shape, HSC Galaxies')
plt.xlim([-3,3])
plt.xlabel('R11')
plt.legend(loc='upper left', fontsize='x-small')
plt.subplot(122)
plt.hist(mcal_gold_matches['R11'][R11_gold_cut*shape_mcal_gold_cut], bins=50, histtype='step', normed=True, label='MCAL Shape Objects')
plt.hist(mcal_hsc_gold_matches['R11'][R11_hsc_cut*star_hsc_cut], bins=50, histtype='step', normed=True, label='All HSC Stars')
plt.hist(mcal_hsc_gold_matches['R11'][R11_hsc_cut*gal_hsc_cut], bins=50, histtype='step', normed=True, label='All HSC Galaxies')
plt.xlim([-3,3])
plt.xlabel('R11')
plt.legend(loc='upper left', fontsize='x-small')
plt.savefig('plots/R11-shapecut-stargal.png', dpi=300, bbox_inches='tight')