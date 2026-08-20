tic; 
for i = 1:1000
    rand(1000,1000)^2;
end
toc;
tic; 
parfor i = 1:1000
    rand(1000,1000)^2;
end
toc;