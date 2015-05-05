type PSF
    psf::ASCIIString
    psfcube::Array{Float64}
    psfavg::Array{Float64}
    nu::Array{Float64}
    nu0::Float64
    mypsf::Array{Float64}
    mypsfadj::Array{Float64}
end

type SKY
    obj::ASCIIString
    sky0::Array{Float64}
    sky::Array{Float64}
    alpha::Array{Float64}
    skyconv::Array{Float64}
    sig::Float64
    var::Float64
    noise::Array{Float64}
    mydata::Array{Float64}
    sumsky2::Array{Float64}
end

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
    s::SharedArray{Float64}
    taus::Array{Float64}
    rhos::Float64
    sh::SharedArray{Float64}

    taup::Array{Float64}
    p::Array{Float64}
    rhop::Float64

    tauv::Array{Float64}
    v::Array{Float64}
    rhov::Float64

    t::Array{Float64}
    taut::Array{Float64}
    rhot::Float64

    wlt::SharedArray{Float64}

    x::SharedArray{Float64}
    Hx::SharedArray{Float64}
    xmm::Array{Float64}

    spectralwlt::Array{Float64}

    fty::Array{Float64}

    μt::Float64
    μv::Float64
    muesp::Float64
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

function init_PSF()
    return PSF("",[],[],[],0.,[],[])
end
function init_SKY()
    return SKY("",[],[],[],[],0.,0.,[],[],[])
end
function init_Algoparam()
    return Algo_param(0,0,0,0,0,0,0)
end
function init_Admmarray()
    return Admm_array([],[],0.,[],[],[],0.,[],[],0.,[],[],0.,[],[],[],[],[],[],0.,0.,0.,0.,0.)
end
function init_TOOLS()
    return TOOLS([],[],[],[],[],[],[],[],[],[])
end
