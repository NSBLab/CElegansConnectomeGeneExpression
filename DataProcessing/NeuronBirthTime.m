% reorder neuron birth time data
birthTimes = GiveMeBirthTimes(C);
birthTimes = celegans279_birth_time(ind);

birthTimeDif = zeros(279,279);
for i=1:279
    for j=1:279
        birthTimeDif(i,j) = abs(birthTimes(j)-birthTimes(i));
    end
end

C.BirthTimeDif = birthTimeDif;
C.BirthTime = birthTimes;

% plot distribution of birth times
figure; histogram(birthTimes, 20); xlabel('birth time (min)'); ylabel('Number of neurons');
% get chemical connectivity matrix
Adj = GiveMeAdj(C,'zeroBinary');
[~,~,deg] = degrees_dir(Adj);
figure; scatter(birthTimes, deg); xlabel('birth time (min)'); ylabel('degree');

% plot the histogram of birth time differences
nanMatrix = tril(nan(279,279));
dif = birthTimeDif+nanMatrix; %take one half of the matrix
figure; histogram(dif); xlabel('birth time difference (min)'); ylabel('Number of neuron pairs');


% plot distribution for neuronal birth time difference for connected links
mask = triu(logical(Adj+Adj')+0);
mask = mask+nanMatrix;
difference = C.BirthTimeDif.*mask;

difference = difference(~isnan(difference));
connection = mask(~isnan(mask));

% combine connection and birth time difference data to select birth time
% differences for existing connections.
data =[difference, connection];
data( ~any(data,2), : ) = [];

figure; histogram(data(:,1), 20); title ('Distribution of neuron birth time differences for connected neurons');
xlabel('Birth time difference (min)'); ylabel('Number of neuron pairs');
