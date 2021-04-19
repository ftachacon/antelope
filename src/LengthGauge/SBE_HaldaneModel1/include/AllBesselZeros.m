%%
clear all
for i=0:4
    zz(i+1,:)=bessel_zeros(1, i, 20001, 1e-8);
    
end

fid = fopen('besselZerosNew.bin', 'w');
fwrite(fid, zz, 'double');
fclose(fid)


%%