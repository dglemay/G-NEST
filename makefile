##########################################################################
# Copyright (C) 2009-2012 William F. Martin and Danielle G. Lemay
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation;
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
##########################################################################

INCLUDE_DIR=`pg_config --includedir`
SRVR_INCLUDE_DIR=`pg_config --includedir-server`
CFLAGS = -g -O3 -std=gnu99 -fpic  -I ${INCLUDE_DIR} -I ${SRVR_INCLUDE_DIR}

TARGETS = gnest_anc  gnest_fill_nh_by_win_size.so  gnest_dot_product.so

%.so : %.o
	ld -shared -o $@ $<

gnest_anc : gnest_anc.o
	gcc $^ -L/usr/lib -fpic -lm -lpq -o $@

all :  ${TARGETS}


clean:
	rm -f ${TARGETS} gnest_anc.o 
