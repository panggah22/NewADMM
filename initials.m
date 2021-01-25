function a = initials(x,steps)

conA = cell(1,steps);
for i = 1:steps
    if i == 1
        conA{i} = x;
    else
        conA{i} = zeros(size(x));      
    end
end
a = cell2mat(conA);

end