B2GAnaFW
========

Analysis framework for Beyond Two Generations (B2G) Physics Analysis Group (PAG) of the Compact Muon Solenoid (CMS) Experiment

Checkout Instructions
=====================

Make a new CMSSW area in 7_3_X

$ cmsrel CMSSW_7_3_0

$ cd CMSSW_7_3_0/src

$ cmsenv

(check out any CMSSW packages needed)

Make an Analysis directory :

$ mkdir Analysis

$ cd Analysis


Clone the github repository

$ git clone https://github.com/cmsb2g/B2GAnaFW.git

$ cd ..

Compile

$ scram b -j 10

Running
=======

The python configuration file for cmsRun is B2GAnaFW/test/b2gedmntuples_cfg.py. It runs on the miniAOD data tier and produces an EDM-ntuple.

The configuration file contians a header explaining usage. Here is an example running on a phys14 top sample:

$ cmsRun b2gedmntuples_cfg.py maxEvts=100 sample="/store/mc/Phys14DR/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/E873348E-BC70-E411-BFA8-0025907B4FD6.root" outputLabel="myoutput"
