####################################
####################################
type PSF
    nu::Array{Float64}
    nu0::Float64
    mypsf::Array{Float64}
    mypsfadj::Array{Float64}
end

type SKY
    sky::Array{Float64}
    alpha::Array{Float64}
    sig::Float64
    var::Float64
    noise::Array{Float64}
    mydata::Array{Float64}
    sumsky2::Array{Float64}
end
####################################
####################################
type PSF_dirty
    nu::Array{Float64}
    nu0::Float64
    mypsf::Array{Float64}
    mypsfadj::Array{Float64}
end

type SKY_dirty
    sky::Array{Float64}
    mydata::Array{Float64}
    sumsky2::Array{Float64}
end
####################################
####################################
type Algo_param
    nspat::Int64
    nfreq::Int64
    nspec::Int64
    nxy::Int64
    niter::Int64
    lastiter::Int64
    nitermax::Int64
end

type Admm_array
    s::Array{Float64}
    taus::Array{Float64}
    rhos::Float64
    sh::Array{Float64}

    taup::Array{Float64}
    p::Array{Float64}
    t::Array{Float64}
    taut::Array{Float64}
    wlt::Array{Float64}
    x::Array{Float64}

    rhop::Float64
    tauv::Array{Float64}
    v::Array{Float64}
    rhov::Float64
    rhot::Float64
    xmm::Array{Float64}

    spectralwlt::Array{Float64}

    fty::Array{Float64}

    μt::Float64
    μv::Float64
    mueps::Float64
    tt::Float64
    mu::Float64
end

type Admm_array_parallel
    s::Array{Float64}
    taus::Array{Float64}
    rhos::Float64
    sh::Array{Float64}

    taup::SharedArray{Float64}
    p::SharedArray{Float64}
    t::SharedArray{Float64}
    taut::SharedArray{Float64}
    wlt::SharedArray{Float64}
    x::SharedArray{Float64}

    rhop::Float64
    tauv::Array{Float64}
    v::Array{Float64}
    rhov::Float64
    rhot::Float64
    xmm::Array{Float64}

    spectralwlt::Array{Float64}

    fty::Array{Float64}

    μt::Float64
    μv::Float64
    mueps::Float64
    tt::Float64
    mu::Float64
end

type TOOLS
    snr::Array{Float64}
    tol1::Array{Float64}
    tol2::Array{Float64}
    tol3::Array{Float64}
    tol4::Array{Float64}
    tol5::Array{Float64}
    err::Array{Float64}
    errorrec::Array{Float64}
    errorest::Array{Float64}
    errorraw::Array{Float64}
end
####################################
####################################
function init_PSF()
    return PSF([],0.,[],[])
end
function init_SKY()
    return SKY([],[],0.,0.,[],[],[])
end
####################################
####################################
function init_PSF_dirty()
    return PSF_dirty([],0.,[],[])
end
function init_SKY_dirty()
    return SKY_dirty([],[],[])
end
####################################
####################################
function init_Algoparam()
    return Algo_param(0,0,0,0,0,0,0)
end
# function init_Admmarray()
#     return Admm_array([],[],0.,[],[],[],0.,[],[],0.,[],[],0.,[],[],[],[],[],0.,0.,0.,0.,0.)
# end

function init_Admmarray(parallel)
    if parallel == "true"
        return Admm_array_parallel([],[],0.,[],[],[],[],[],[],[],0.,[],[],0.,0.,[],[],[],0.,0.,0.,0.,0.)
    else
        return Admm_array([],[],0.,[],[],[],[],[],[],[],0.,[],[],0.,0.,[],[],[],0.,0.,0.,0.,0.)
    end
end

function init_TOOLS()
    return TOOLS([],[],[],[],[],[],[],[],[],[])
end
