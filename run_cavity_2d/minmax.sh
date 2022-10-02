#!/usr/bin/env bash
more nohup.out | grep U,V,W | grep min | cut -c12-43 > uvmin
more nohup.out | grep U,V,W | grep max | cut -c12-43 > uvmax
more nohup.out | grep Phi1 | cut -c15- > phiminmax
