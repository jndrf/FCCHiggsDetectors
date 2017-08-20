#!/bin/sh

for f in *.py; do
    name=$(basename $f .py)
    sed "s|from heppy.papas.detectors.CMS import CMS|from heppy.papas.detectors.FCCHiggsDetectors.$name import CMS|" <$HEPPY/test/papas_cfg.py > config/cfg_$f
done