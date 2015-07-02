#!/bin/bash
IN=$1
OUT=$2
julia /home/mirco_local/glob/GaussDCA.jl/rungDCA.jl $IN $OUT > ${IN%.trimmed}.gneff
