function matmerge(dest,file1,file2)
%MATMERGE Merge .mat-files.
%   MATMERGE(DEST,FILE1,FILE2) Loads the two files FILE1 and FILE2 and
%   examines the variables. The intersection of variable names is
%   determined and the values of all those variables are concatenated
%   along the first non-matching dimension (or last dimension). The
%   resulting variables are then stored in the created file DEST. It
%   is considered an error if DEST exists.
%
%   MATMERGE(DEST,FILE1) allows for an existing file DEST. Hence this
%   syntax can be thought of as DEST += FILE1.
%
%   MATMERGE is a simple design to aid with distributed computations. No
%   attempt in merging non-array data or making the data unique has been
%   made.

% S. Engblom 2012-04-01

% new file or appending to dest
if ~strcmp(dest,'.mat')
  dest = [dest '.mat'];
end
if nargin > 2
  if exist(dest,'file')
    error('matmerge:e1','Will not overwrite existing file.');
  end
else
  if ~exist(dest,'file')
    error('matmerge:e2','Cannot merge with non-existing file.');
  end
  file2 = file1;
  file1 = dest;
end

% load files
s1 = load(file1);
n1 = fieldnames(s1);
s2 = load(file2);
n2 = fieldnames(s2);

% intersection of variable names
n3 = fsetop('intersect',n1,n2);
s3 = struct;

% merge data
for i = 1:size(n3,2)
  temp1 = s1.(n3{i});
  temp2 = s2.(n3{i});
  effdim = max(tndims(temp1),tndims(temp2));
  sz1 = tsize(temp1,1:effdim);
  sz2 = tsize(temp2,1:effdim);
  catdim = find([sz1 0] ~= [sz2 1],1,'first');
  % this wont work if temp1 and temp2 doesn't match
  s3.(n3{i}) = cat(catdim,temp1,temp2);
end

save(dest,'-struct','s3');
