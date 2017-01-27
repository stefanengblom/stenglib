function hash = perfecthash(s)
%PERFECTHASH Perfect hash-function from strings.
%   HASH = PERFECTHASH(S) with S a cell-vector of strings tries to 
%   determine a perfect hash-function. A few simple models for such 
%   functions are tried and any function determined in this way is 
%   returned as a string HASH in C pseudo-code.
%
%   Note: PERFECTHASH only indexes characters with index <= maximum length 
%   of any of the strings. You must explicitly increase the length of the 
%   strings to index other characters. See the second example below.
%
%   Example:
%      s = {'check','intersect','setdiff','setxor', ...
%           'union','unique','ismember'};
%      perfecthash(s)
%
%      % case when increasing the length of the strings might help
%      s = {'if' 'else' 'elseif' 'while' 'for' 'switch' 'case' ...
%           'end' 'return' 'break' 'continue' 'function'};
%      perfecthash(s) % void
%
%      % add strlen(s{i}) and ending '\0'
%      sz = cellfun('prodofsize',s);
%      for i = 1:numel(s), s{i} = [char(sz(i)) s{i} char(0)]; end
%      perfecthash(s)
%      % (here str[0] is strlen(str))

% S. Engblom 2016-10-07 (Improved handling of non-printable characters)
% S. Engblom 2011-06-27 (enum+char*[] output, bug for 3-ary hashes fixed)
% S. Engblom 2008-04-25 (Minor revision)
% S. Engblom 2007-05-25

% to be built
hash = '';

% construct character array c with strings as rows
c = char(s);
slen = min(cellfun('prodofsize',s));
nums = size(c,1);
% consistently using doubles (not optimal -- but convenient)
c = double(c(:,1:slen));

% enums are formatted in printable uppercases, strings are case-sensitives
S = s;
for i = 1:nums
  prints1 = find(33 <= s{i} & s{i} <= 126);
  prints2 = ('a' <= s{i} & s{i} <= 'z') | ...
            ('A' <= s{i} & s{i} <= 'Z') | ...
            ('0' <= s{i} & s{i} <= '9');
  s{i} = s{i}(prints1); % be more generous here
  S{i} = upper(S{i});
  S{i}(~prints2) = '_';
end

% try these sizes for the hash-tables
siz = nums:2*nums;

% factors to try
muls = 1:20;

% 1st order hash-functions
for sz = siz
  cc = mod(c,sz);
  ii = all(fsparse(cc+1,1:slen,1) <= 1,1);
  if any(ii)
    for i = find(ii)
      % output function
      str = ['typedef enum {' frepmat('%s = %d,',[1 size(c,1)])];
      str = [str(1:end-1) '} PROP;\n'];
      [foo,is] = sort(cc(:,i));
      t = [S(is); mat2cell(cc(is,i)',1,ones(1,size(c,1)))];
      hash = [hash sprintf(str,t{:})];

      t = frepmat({'NULL%s,'},[1 sz]);
      t(cc(:,i)+1) = {'&get%s,'};
      t = cell2mat(t);
      str = ['static const getfun_t getval[] = {' t];
      str = [str(1:end-1) '};\n'];
      t = frepmat({''},[1 sz]);
      t(cc(:,i)+1) = S;
      hash = [hash sprintf(str,t{:})];

      t = frepmat({'NULL%s,'},[1 sz]);
      t(cc(:,i)+1) = {'&%s,'};
      t = cell2mat(t);
      str = ['void *val[] = {' t];
      str = [str(1:end-1) '};\n'];
      t = frepmat({''},[1 sz]);
      t(cc(:,i)+1) = S;
      hash = [hash sprintf(str,t{:})];

      str = ['static const char *tab[] = {' frepmat('"%s",',[1 sz])];
      str = [str(1:end-1) '};\n'];
      t = frepmat({''},[1 sz]);
      t(cc(:,i)+1) = s;
      hash = [hash sprintf(str,t{:})];
      hash = [hash sprintf('hash(str) = str[%d] %% %d;\n\n',i-1,sz)];
    end
    % no point in continuing here
    siz = nums:sz-1;
    break;
  end
end

% 2nd order functions

% *** use of symmetricity could reduce time by a factor of 2
found = false;
for sz = siz
  for mul1 = muls
    for mul2 = muls
      cc = reshape(tsum(mul1*c,mul2*c,[1 2],[1 3]),nums,[]);
      cc = mod(cc,sz);
      ii = all(fsparse(cc+1,1:slen^2,1) <= 1,1);
      if any(ii)
	[i1,i2] = ind2sub([slen slen],find(ii));
	ix = i1 <= i2;
	if any(ix)
	  i1 = i1(ix); i2 = i2(ix);
          for i = 1:numel(i1)
            j = (i2(i)-1)*slen+i1(i);
            str = ['typedef enum {' frepmat('%s = %d,',[1 size(c,1)])];
            str = [str(1:end-1) '} PROP;\n'];
            [foo,is] = sort(cc(:,j));
            t = [S(is); ...
                 mat2cell(cc(is,j)',1,ones(1,size(c,1)))];
            hash = [hash sprintf(str,t{:})];
            
            t = frepmat({'NULL%s,'},[1 sz]);
            t(cc(:,j)+1) = {'&get%s,'};
            t = cell2mat(t);
            str = ['static const getfun_t getval[] = {' t];
            str = [str(1:end-1) '};\n'];
            t = frepmat({''},[1 sz]);
            t(cc(:,j)+1) = S;
            hash = [hash sprintf(str,t{:})];

            t = frepmat({'NULL%s,'},[1 sz]);
            t(cc(:,j)+1) = {'&%s,'};
            t = cell2mat(t);
            str = ['void *val[] = {' t];
            str = [str(1:end-1) '};\n'];
            t = frepmat({''},[1 sz]);
            t(cc(:,j)+1) = S;
            hash = [hash sprintf(str,t{:})];

            str = ['static const char *tab[] = {' frepmat('"%s",',[1 sz])];
            str = [str(1:end-1) '};\n'];
            t = frepmat({''},[1 sz]);
            t(cc(:,j)+1) = s;
            hash = [hash sprintf(str,t{:})];
            hash = [hash ...
                     sprintf(['hash(str) = ' ...
		       '((str[%d] * %d)+(str[%d] * %d)) %% %d;\n\n'], ...
	               i1(i)-1,mul1,i2(i)-1,mul2,sz)];
          end
          siz = nums:sz-1;
          found = true;
          break;
        end
      end
    end
    if found, break; end
  end
  if found, break; end
end

% 3rd order functions

% *** use of symmetricity could reduce time by a factor of 6
for sz = siz
  for mul1 = muls
    for mul2 = muls
      for mul3 = muls
        cc = reshape(tsum(mul1*c,mul2*c,[1 2],[1 3]),nums,[]);
        cc = reshape(tsum(cc,mul3*c,[1 2],[1 3]),nums,[]);
        cc = mod(cc,sz);
        ii = all(fsparse(cc+1,1:slen^3,1) <= 1,1);
        if any(ii)
	  [i1,i2,i3] = ind2sub([slen slen slen],find(ii));
          ix = i1 <= i2 & i2 <= i3;
          if any(ix)
            i1 = i1(ix); i2 = i2(ix); i3 = i3(ix);
            for i = 1:numel(i1)
              j = ((i3(i)-1)*slen+(i2(i)-1))*slen+i1(i);
              str = ['typedef enum {' frepmat('%s = %d,',[1 size(c,1)])];
              str = [str(1:end-1) '} PROP;\n'];
              [foo,is] = sort(cc(:,j));
              t = [S(is); mat2cell(cc(is,j)',1,ones(1,size(c,1)))];
              hash = [hash sprintf(str,t{:})];

              t = frepmat({'NULL%s,'},[1 sz]);
              t(cc(:,j)+1) = {'&get%s,'};
              t = cell2mat(t);
              str = ['static const getfun_t getval[] = {' t];
              str = [str(1:end-1) '};\n'];
              t = frepmat({''},[1 sz]);
              t(cc(:,j)+1) = S;
              hash = [hash sprintf(str,t{:})];

              t = frepmat({'NULL%s,'},[1 sz]);
              t(cc(:,j)+1) = {'&%s,'};
              t = cell2mat(t);
              str = ['void *val[] = {' t];
              str = [str(1:end-1) '};\n'];
              t = frepmat({''},[1 sz]);
              t(cc(:,j)+1) = S;
              hash = [hash sprintf(str,t{:})];

              str = ['static const char *tab[] = {' frepmat('"%s",',[1 sz])];
              str = [str(1:end-1) '};\n'];
              t = frepmat({''},[1 sz]);
              t(cc(:,j)+1) = s;
              hash = [hash sprintf(str,t{:})];

              hash = [hash ...
                       sprintf(['hash(str) = ' ...
		         '((str[%d] * %d)+(str[%d] * %d)+(str[%d] * %d)) ' ...
                         '%% %d;\n'], ...
		         i1(i)-1,mul1,i2(i)-1,mul2,i3(i)-1,mul3,sz)];
            end
            return;
          end
        end
      end
    end
  end
end
