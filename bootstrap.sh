#!/bin/sh -e

aclocal
autoconf
autoheader
automake --add-missing
