#!/bin/bash
pushd "$ROOT/scripts/02_rmsf"

./rmsf_test.r
./rmsf_ratios.r

popd