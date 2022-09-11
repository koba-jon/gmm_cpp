#!/bin/bash

DATA='toy'

./GMM \
    --dataset ${DATA} \
    --D 3 \
    --K 8
