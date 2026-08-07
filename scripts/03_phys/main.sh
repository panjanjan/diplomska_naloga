#!/bin/bash
pushd "$ROOT/scripts/03_phys"

./com_distances_bio3d.r
./principal_axes_of_inertia.r
./psd.r

popd