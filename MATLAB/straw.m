function varargout = straw(norm, file, chr1, chr2, unit, resolution, matrixType)
%STRAW Query a .hic file through the native libstraw MEX binding.
if nargin < 7
    matrixType = 'observed';
end
[varargout{1:nargout}] = straw_mex(norm, file, chr1, chr2, unit, resolution, matrixType);
end
