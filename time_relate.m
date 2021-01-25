function a = time_relate(x_before,x_after,steps)

con = cell(steps-1,steps);

for i = 1:steps-1
    con{i,i} = x_before;
    con{i,i+1} = x_after;
end
for i = 1:steps-1
    for j = 1:steps
        if isempty(con{i,j})
            con{i,j} = zeros(size(x_before));
        end
    end
end
a = cell2mat(con);

end