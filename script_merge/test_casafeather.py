
importfits('HC3N_natural_cubek_65kms.fits', 'HC3N_natural_cubek_65kms.image', overwrite=True)
importfits('HC3N_natural_tpcube_k_spatialandspectralregrid_65kms.fits', 'HC3N_natural_tpcube_k_spatialandspectralregrid_65kms.image', overwrite=True)
importfits('HC3N_natural_cubeJy_65kms.fits', 'HC3N_natural_cubeJy_65kms.image', overwrite=True)
importfits('HC3N_natural_tpcube_Jy_rg_65kms.fits', 'HC3N_natural_tpcube_Jy_rg_65kms.image', overwrite=True)


print("hc3n_feather_comb_test_kelvin_default")
os.system('rm -rf hc3n_feather_comb_test_kelvin_default')
feather(imagename='hc3n_feather_comb_test_kelvin_default',
        highres='HC3N_natural_cubek_65kms.image',
        lowres='HC3N_natural_tpcube_k_spatialandspectralregrid_65kms.image',
        )
exportfits('hc3n_feather_comb_test_kelvin_default',
           'hc3n_feather_comb_test_kelvin_default.fits', overwrite=True)

print("hc3n_feather_comb_test_kelvin_sd10")
os.system('rm -rf hc3n_feather_comb_test_kelvin_sd10')
feather(imagename='hc3n_feather_comb_test_kelvin_sd10',
        highres='HC3N_natural_cubek_65kms.image',
        lowres='HC3N_natural_tpcube_k_spatialandspectralregrid_65kms.image',
        sdfactor=5591,
        )
exportfits('hc3n_feather_comb_test_kelvin_sd10',
           'hc3n_feather_comb_test_kelvin_sd10.fits', overwrite=True)

print("hc3n_feather_comb_test_jy_default")
os.system('rm -rf hc3n_feather_comb_test_jy_default')
feather(imagename='hc3n_feather_comb_test_jy_default',
        highres='HC3N_natural_cubeJy_65kms.image',
        lowres='HC3N_natural_tpcube_Jy_rg_65kms.image',
        )
exportfits('hc3n_feather_comb_test_jy_default',
           'hc3n_feather_comb_test_jy_default.fits', overwrite=True)
