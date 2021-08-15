k=9
function test(x,y)
    x+y
end

function trysum(i)
    sum(j->test(i,j),1:(k+1))
end
