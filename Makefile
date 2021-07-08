ALL: make_wssd_sunk make_ssf_DTS_caller

make_wssd_sunk:
	pushd wssd_sunk; make; popd

make_ssf_DTS_caller:
	pushd ssf_DTS_caller; make; popd
