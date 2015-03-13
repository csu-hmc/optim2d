function JRF = getJRF(R_JRF, L_JRF)

JRF = zeros(size(R_JRF));
for i = 1:size(R_JRF,2)
    JRF(1,i) = norm(R_JRF(5:6,i));
    JRF(2,i) = norm(R_JRF(3:4,i));
    JRF(3,i) = norm(R_JRF(1:2,i));
    JRF(4,i) = norm(L_JRF(5:6,i));
    JRF(5,i) = norm(L_JRF(3:4,i));
    JRF(6,i) = norm(L_JRF(1:2,i));
end