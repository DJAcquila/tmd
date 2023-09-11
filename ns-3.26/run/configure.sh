#!/bin/bash

./waf clean
./waf -d debug --enable-examples --enable-tests configure
