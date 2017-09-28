function birthTimes = GiveMeBirthTimes(C)
% Loads birth time data and matches to neurons
%-------------------------------------------------------------------------------

load('celegans279_BT');
Names = C.RegionAcronyms;
[together,ind] = intersect(celegans279labels,Names);
if length(together) < length(Names)
    error('Error matching birth times to neurons');
end
% Labels = celegans279labels(ind);
birthTimes = celegans279_birth_time(ind);

end
