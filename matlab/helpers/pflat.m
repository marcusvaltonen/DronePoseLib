function out=pflat(x)
out = x./repmat(x(end,:),[size(x,1) 1]);
end
