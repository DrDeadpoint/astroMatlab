function [elem] = osculating_elements(traj_in,primary)
% [elem] = osculating_elements(traj_in,primary)
%
% RJ has a computationally efficient method
% this function calculates osculating keplerian elements and puts them into a structure

traj_inert = traj_in.changeFrame([primary 'centinert']);
trajLength = traj_inert.getLength;
elem = eci_elem(traj_inert); %contains first index of states
sma = zeros(1,trajLength);
ecc = zeros(1,trajLength);
inc = zeros(1,trajLength);
raan = zeros(1,trajLength);
argP = zeros(1,trajLength);
trueAnom = zeros(1,trajLength);
parfor i = 1:trajLength
    traj = traj_inert.extractSubTraj(i,trajLength);
    thisElem = eci_elem(traj);
    sma(i) = thisElem.sma.value;
    ecc(i) = thisElem.ecc.value;
    inc(i) = thisElem.inc.value;
    raan(i) = thisElem.raan.value;
    argP(i) = thisElem.argP.value;
    trueAnom(i) = thisElem.trueAnom.value;
end
elem.sma.value = sma;
elem.ecc.value = ecc;
elem.inc.value = inc;
elem.raan.value = raan;
elem.argP.value = argP;
elem.trueAnom.value = trueAnom;
end