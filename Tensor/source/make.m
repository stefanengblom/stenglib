function make
%MAKE Makefile for TENSOR.

% S. Engblom 2019-01-23 (mexmaci64, mexa64, 9.6)
% S. Engblom 2015-03-20 (mexa64, 8.4)
% S. Engblom 2015-01-19 (mexmaci64, 8.4)
% S. Engblom 2012-04-16 (mexmaci64, 7.11)
% S. Engblom 2011-04-17 (mexmaci64, 7.10)
% S. Engblom 2011-03-07 (mexa64, 7.11)
% S. Engblom 2010-09-23 (mexs64, 7.7)
% S. Engblom 2010-02-02 (mexa64, 7.8)
% S. Engblom 2010-01-12 (mexmaci)
% S. Engblom 2007-05-17 (mexs64)
% S. Engblom 2006-11-09 (mexa64)
% S. Engblom 2005-04-10 (mexmac)

% Use '-DBLASINT=size_t' for the (bad!) platforms where the 'int' in
% the declaration of BLAS subroutines is in fact a 'size_t' and
% sizeof(size_t) > sizeof(int).

s = pwd;
mx = mexext;
ver = version;

if strcmp(mx,'mexglx')
  if ~strncmp(version,'7.5',3) && ~strncmp(version,'7.8',3)
    warning(['Extension .' mexext [' tested with Matlab version(s) ' ...
                        '7.5 and 7.8 only.']]);
  end
  mex(['CFLAGS=-fPIC -fno-omit-frame-pointer -std=c99 ' ...
       '-D_GNU_SOURCE -pthread -fexceptions'], ...
      '-outdir',s,[s '/source/tndims.c']);
  mex(['CFLAGS=-fPIC -fno-omit-frame-pointer -std=c99 ' ...
       '-D_GNU_SOURCE -pthread -fexceptions'], ...
      '-outdir',s,[s '/source/tsize.c']);
  mex(['CFLAGS=-fPIC -fno-omit-frame-pointer -std=c99 ' ...
       '-D_GNU_SOURCE -pthread -fexceptions'], ...
      '-outdir',s,[s '/source/tsum.c']);
  mex(['CFLAGS=-fPIC -fno-omit-frame-pointer -std=c99 ' ...
       '-D_GNU_SOURCE -pthread -fexceptions'], ...
      '-outdir',s,'-lmwblas',[s '/source/tprod.c']);
elseif strcmp(mx,'mexa64')
  v = version;
  if v(1) == '7'
    if ~strncmp(version,'7.2',3) && ~strncmp(version,'7.8',3) && ...
          ~strncmp(version,'7.11',4) && ~strncmp(version,'7.13',4)
      warning(['Extension .' mexext [' tested with Matlab version(s) ' ...
                          '7.2, 7.8, 7.11 and 7.13 only.']]);
    end
    if ~strncmp(version,'7.11',4)
      mex(['CFLAGS=-fPIC -fno-omit-frame-pointer -std=c99 ' ...
           '-D_GNU_SOURCE -pthread -fexceptions'], ...
          '-outdir',s,[s '/source/tndims.c']);
      mex(['CFLAGS=-fPIC -fno-omit-frame-pointer -std=c99 ' ...
           '-D_GNU_SOURCE -pthread -fexceptions'], ...
          '-outdir',s,[s '/source/tsize.c']);
      mex(['CFLAGS=-fPIC -fno-omit-frame-pointer -std=c99 ' ...
           '-D_GNU_SOURCE -pthread -fexceptions'], ...
          '-outdir',s,[s '/source/tsum.c']);
      if strncmp(version,'7.2',3)
        mex(['CFLAGS=-fPIC -fno-omit-frame-pointer -std=c99 ' ...
             '-D_GNU_SOURCE -pthread -fexceptions'], ...
            '-outdir',s,[s '/source/tprod.c']);
      else
        mex(['CFLAGS=-fPIC -fno-omit-frame-pointer -std=c99 ' ...
             '-D_GNU_SOURCE -pthread -fexceptions -DBLASINT=size_t'], ...
            '-outdir',s,'-lmwblas',[s '/source/tprod.c']);
      end
    else
      % apparently, the linker path is not properly set up on 7.11:
      mex(['CFLAGS=-fPIC -fno-omit-frame-pointer -std=c99 -O3 ' ...
           '-D_GNU_SOURCE -pthread -fexceptions'], ...
          ['-L' matlabroot '/sys/os/glnxa64'], ...
          '-outdir',s,[s '/source/tndims.c']);
      mex(['CFLAGS=-fPIC -fno-omit-frame-pointer -std=c99 -O3 ' ...
           '-D_GNU_SOURCE -pthread -fexceptions'], ...
          ['-L' matlabroot '/sys/os/glnxa64'], ...
          '-outdir',s,[s '/source/tsize.c']);
      mex(['CFLAGS=-fPIC -fno-omit-frame-pointer -std=c99 -O3 ' ...
           '-D_GNU_SOURCE -pthread -fexceptions'], ...
          ['-L' matlabroot '/sys/os/glnxa64'], ...
          '-outdir',s,[s '/source/tsum.c']);
      mex(['CFLAGS=-fPIC -fno-omit-frame-pointer -std=c99 -O3 ' ...
           '-D_GNU_SOURCE -pthread -fexceptions -DBLASINT=size_t'], ...
          ['-L' matlabroot '/sys/os/glnxa64'], ...
          '-outdir',s,'-lmwblas',[s '/source/tprod.c']);
    end
  else
    if ~strncmp(version,'8.4',3) && ~strncmp(version,'9.6',3)
      warning(['Extension .' mexext ' tested with Matlab version(s) ' ...
               '8.4 and 9.6 only.']);
    end
      
    % apparently, the linker path is not properly set up on 8.4 (also a
    % soft link libstdc++.so inside [matlabroot '/sys/os/glnxa64'] is
    % required to point to the correct shared library, in this case
    % libstdc++.so.6.0.17)
    mex('CFLAGS=-fPIC -std=c99 -O3',['-L' matlabroot '/sys/os/glnxa64'], ...
        '-outdir',s,[s '/source/tndims.c']);
    mex('CFLAGS=-fPIC -std=c99 -O3',['-L' matlabroot '/sys/os/glnxa64'], ...
        '-outdir',s,[s '/source/tsize.c']);
    mex('CFLAGS=-fPIC -std=c99 -O3',['-L' matlabroot '/sys/os/glnxa64'], ...
        '-outdir',s,[s '/source/tsum.c']);
    mex('CFLAGS=-fPIC -std=c99 -O3 -DBLASINT=size_t', ...
        ['-L' matlabroot '/sys/os/glnxa64'], ...
        '-outdir',s,'-lmwblas',[s '/source/tprod.c']);
  end
elseif strcmp(mx,'mexmac')
  if ~strncmp(version,'7.0',3)
    warning(['Extension .' mexext ' tested with Matlab version(s) 7.0 only.']);
  end
  mex('CC=gcc -std=c99','-outdir',s,[s '/source/tndims.c']);
  mex('CC=gcc -std=c99','-outdir',s,[s '/source/tsize.c']);
  mex('CC=gcc -std=c99','-outdir',s,[s '/source/tsum.c']);
  mex('CC=gcc -std=c99','-outdir',s,[s '/source/tprod.c']);
elseif strcmp(mx,'mexmaci')
  if ~strncmp(version,'7.8',3)
    warning(['Extension .' mexext ' tested with Matlab version(s) 7.8 only.']);
  end
  mex('CC=gcc -std=c99 -fast','-outdir',s,[s '/source/tndims.c']);
  mex('CC=gcc -std=c99 -fast','-outdir',s,[s '/source/tsize.c']);
  mex('CC=gcc -std=c99 -fast','-outdir',s,[s '/source/tsum.c']);
  mex('CC=gcc -std=c99 -fast','-outdir',s,'-lmwblas', ...
      [s '/source/tprod.c']);
elseif strcmp(mx,'mexmaci64')
  v = version;
  if v(1) == '7'
    if ~strncmp(version,'7.10',4) && ~strncmp(version,'7.11',4) && ...
          ~strncmp(version,'7.14',4)
      warning(['Extension .' mexext ' tested with Matlab version(s) ' ...
               '7.10 and 7.11 only.']);
    end
    mex('CC=gcc -std=c99 -fast','-outdir',s,[s '/source/tndims.c']);
    mex('CC=gcc -std=c99 -fast','-outdir',s,[s '/source/tsize.c']);
    mex('CC=gcc -std=c99 -fast','-outdir',s,[s '/source/tsum.c']);
    mex('CC=gcc -std=c99 -fast -DBLASINT=size_t', ...
        '-outdir',s,'-lmwblas', ...
        [s '/source/tprod.c']);
  else
    if ~strncmp(version,'8.4',3) && ~strncmp(version,'9.6',3)
      warning(['Extension .' mexext ' tested with Matlab version(s) ' ...
               '8.4 and 9.6 only.']);
    end
    mex('CFLAGS= -std=c99','-outdir',s,[s '/source/tndims.c']);
    mex('CFLAGS=-Wno-logical-op-parentheses -std=c99','-outdir',s,[s '/source/tsize.c']);
    mex('CFLAGS=-Wno-logical-op-parentheses -std=c99','-outdir',s,[s '/source/tsum.c']);
    mex('CFLAGS= -std=c99 -DBLASINT=size_t', ...
        '-outdir',s,'-lmwblas', ...
        [s '/source/tprod.c']);
  end
elseif strcmp(mx,'mexs64')
  if ~strncmp(version,'7.7',3)
    warning(['Extension .' mexext ' tested with Matlab version(s) 7.7 only.']);
  end
  mex('-outdir',s,[s '/source/tndims.c']);
  mex('-outdir',s,[s '/source/tsize.c']);
  mex('-outdir',s,[s '/source/tsum.c']);
  mex('-outdir',s,[s '/source/tprod.c']); 
else
  warning('New platform. Trying default make.');
  mex('-outdir',s,[s '/source/tndims.c']);
  mex('-outdir',s,[s '/source/tsize.c']);
  mex('-outdir',s,[s '/source/tsum.c']);
  mex('-outdir',s,[s '/source/tprod.c']); 
end
