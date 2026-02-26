#!/bin/bash

set -ex

export CAPYRXSPACE="$PWD"
export WORKSPACE="$PWD/../workspace"
mkdir -p "$WORKSPACE"
cd "$WORKSPACE"

# Check out Cactus
wget https://raw.githubusercontent.com/gridaphobe/CRL/master/GetComponents
chmod a+x GetComponents
./GetComponents --no-parallel --shallow "$CAPYRXSPACE/scripts/capyrx.th"

cd Cactus

# Create a link to the CapyrX repository
ln -s "$CAPYRXSPACE" repos
# Create links for the CapyrX thorns
mkdir -p arrangements/CapyrX
pushd arrangements/CapyrX
ln -s ../../repos/CapyrX/* .
popd
