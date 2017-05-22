import runpy
result = runpy.run_path('core_photometry.py', run_name='__main__')
assert 'tbl' in result
runpy.run_path('merge_core_masers.py', run_name='__main__')
runpy.run_path('main_photometry_table.py', run_name='__main__')
runpy.run_path('stellar_mass_estimates.py', run_name='__main__')
