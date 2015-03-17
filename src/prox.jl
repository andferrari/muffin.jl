# proximity operators
function prox_visi(y_est::Array)
    return (2.*visi + rho_y.*y_est)./(rho_y + 2)
end

function prox_u(u::Array,μ::const)
    return (max(1-μ./abs(u),0).*u)
end
