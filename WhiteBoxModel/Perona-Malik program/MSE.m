function Error = MSE(imageA, imageB)
  if(length(size(imageA)) != 2 && length(size(imageB)) != 2)
    break;
  end
  Error = sum(sum((imageA - imageB) .^2))
end