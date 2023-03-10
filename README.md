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
programs:
```
./generateSideband savesb_[era]_b[# of bin].root
it creates 
- a directory  GenBin[# of bin]_[era] 
- generates 100 toys using the sb model found in  savesb_[era]_b[# of bin].root and it save them in the direcory GenBin[# of bin]_[era]
- save histograms and canvas in this file: generate-savesb_[era]_b[# of bin].root
```
```
./readSideband savesb_[era]_b[# of bin].root
- read the sb model and produce a plot with the angular projections of the  model 
```
