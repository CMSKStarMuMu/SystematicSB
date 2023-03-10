# SystematicSB
## Setup working area
```sh
export SCRAM_ARCH=slc7_amd64_gcc820
cmsrel CMSSW_10_4_0
cd CMSSW_10_4_0/src/ && cmsenv && cd ../..
```
Clone this branch in the working directory:
```sh
git clone -b master git@github.com:CMSKStarMuMu/SystematicSB.git
cd SystematicSB
```
to compile
```
make
```
```
generateSideband
```
