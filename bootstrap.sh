#!/bin/sh -e

aclocal -I m4 --install
autoconf
autoheader
automake --add-missing
