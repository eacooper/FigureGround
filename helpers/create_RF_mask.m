function RFmask = create_RF_mask(r,c,RF_size_pix,col_vals,row_vals)
%
% take a pixel location and RF distanceter and the indices of the rows/columns in a matrix and
% create a binary mask for pixels within the RF

% initialize mask
RFmask = zeros(size(col_vals));

% get distance of all pixels from the RF center
dist = sqrt((r - row_vals).^2 + (c - col_vals).^2);

% find pixels within RF
RFmask(dist <= round(RF_size_pix/2)) = 1;

