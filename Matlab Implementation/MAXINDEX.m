function[max_index] = MAXINDEX(j,arr)

shape = size(arr);
max_index = [1,1];

for i=1:shape(2)
    for k=1:shape(3)
        if (arr(j, max_index(1),max_index(2))< arr(j,i,k))
            max_index = [i, k];
        end
    end
end
