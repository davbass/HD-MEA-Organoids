"""
    A collection of functions for analysis of BRW (HDF5) files generated with the BrainWave program from the
    company 3Brain.
    Laboratory 19 of the CINVESTAV in charge of Dr. Rafael Gutierrez Aguilar.
    Work developed by Isabel Romero-Maldonado (2020 - )
    isabelrm.biofisica@gmail.com
    https://github.com/LBitn
"""

module AccuraTools

using HDF5
using JLD
using Dates
using DelimitedFiles
using BinningAnalysis
using Distributions
using HistogramThresholding
using ImageContrastAdjustment
using Suppressor
using Plots
using Statistics
using DSP
using StatsBase

export VariablesBRW4, Get_Chunks, div_ab, searchdir, VoltageConvertMeanΔxCI, StaticThresholdEvents,
       Thresholding, ThresholdingPlots, Z0, ZW, Zplot, STExbin, STExChannel, FiltMUAremez, FilteringMatrix,
       ms2frames, ΔV, VoltageConvert, PositiveSaturation, donoho, FillingHolesCrux, Get_Groups, neighborgs,
       FigureGroups, Cut_Spikes, FiltM, PTXchannel, PosT;

# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
    VariablesBRW4( FILEBRW::String, P = false ) -> Variables::Dict{Any, Any}
        Extracts all the metadata in the BRW file. Saves the variables into the parent folder of the file.
        If P = true, also deletes the entries of the dictionary that appear to be null or empty
        using HDF5, JLD, DelimitedFiles, AccuraTools.ExperimentSettings2Dict, AccuraTools.ExperimentDate2String
"""
function VariablesBRW4( FILEBRW::String, P = false )
    BRW = h5open( FILEBRW, "r" );
    Raw = [ ]
    GROUPRAW = [ ]
    for i in BRW
        if haskey( i, "Raw" )
            Raw = string( i );
            GROUPRAW = split( Raw, " " )[ 2 ];
            Raw = string( GROUPRAW, "/Raw" );
        end
    end
    Groups = keys( BRW );
    Attributes = attributes( BRW );
    GROUPRAW = replace( GROUPRAW, "/" => "" );
    D = Dict( );
    for g in Groups
        aux = attributes( BRW[ g ] );
        for j = 1:length( aux );
            D[ string( g, "/", keys( aux )[ j ] ) ] = read( aux[ keys( aux )[ j ] ] );
        end
    end
    for i in keys( Attributes )
        D[ i ] = read( open_attribute( BRW, i ) );
    end
    D[ "Raw" ] = Raw;
    aux = keys( BRW[ GROUPRAW ] )[ keys( BRW[ GROUPRAW ] ) .!= "Raw" ];
    for i in aux
        D[ i ] = read( BRW[ string( GROUPRAW, "/", i ) ] )
    end
    Groups_NORaw = Groups[ Groups .!= GROUPRAW ];
    for g in Groups_NORaw
        D[ g ] = read( BRW[ g ] );
    end
    BRWsize = ( ( stat( FILEBRW ).size ) / 1000000 ) / 1024;
    D[ "BRWname" ] = BRW.filename;
    D[ "BRWsize"] = BRWsize;
    D = ExperimentSettings2Dict( D );
    D = ExperimentDate2String( D );
    D[ "nChs" ] = length( D[ "StoredChIdxs" ] );
    if P == true
        X = String.( keys( D ) )[ values( D ) .== "null" ];
        for i in X
            delete!( D, i );
        end
        X = [ ];
        for i in keys( D )
            try
                if isempty( D[ i ] )
                    push!( X, i );
                end
            catch e
            end
        end
        for i in X
            delete!( D, i );
        end
    end
    x = [ ];
    for i in values( D )
    if length( i ) <= 50; push!( x, 1 ); else; push!( x, 0 ); end
    end
    NK = string.( keys( D ) )[ Bool.( x ) ];
    TXT = Dict( ); for i in NK; TXT[ i ] = D[ i ]; end
    PATHMAIN = split( FILEBRW, "." )[ 1 ]; mkpath( PATHMAIN );
    PATHINFO = joinpath( PATHMAIN, "Info" ); mkpath( PATHINFO );
    FILEVARS = joinpath( PATHINFO, "Variables.jld" );
    FILEVARSTXT = joinpath( PATHINFO, "Variables.txt" );

    writedlm( FILEVARSTXT, TXT );
    save( FILEVARS, "Variables", D );

    println( "You are now working on the new main path: ", PATHMAIN );
    println( basename( BRW.filename ), " : ", D[ "Description" ] );
    println( "HDF5 file size: $BRWsize GB" );
    cd( PATHMAIN )

    close( BRW )
    return D

end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
    ExperimentSettings2Dict( Variables::Dict{Any, Any} ) -> Variables::Dict{Any, Any}
        Extracting the contents of the ExperimentSettings dictionary
"""
function ExperimentSettings2Dict( Variables::Dict{Any, Any} )
    ExperimentSettings = Variables[ "ExperimentSettings" ][ 1 ];
    t = split( ExperimentSettings,"\r\n" );
    t = replace.( t, "  " => "", "{" => "", "}" => "", '"' => "" );
    x = [ ];
    for i = 1:length( t )
        if !isempty( collect( eachmatch( r"[a-z]", t[ i ] ) ) )
            push!( x, true );
        else
            push!( x, false );
        end
    end
    t = t[ Bool.( x ) ]; t = split.( t, ": " );
    D = Dict( );
    for i in t
        if !( i[ 2 ] == "" )
            aux = i[ 2 ];
            try
                aux = replace( i[ 2 ], "," => " ", "[" => "", "[]" => "","]" => "", " " => "" )
            catch e
            end
            if ( aux != "" ) && ( aux != " " )
                aux = aux;
                try
                    aux = parse( Float64, aux );
                catch
                    aux = replace( aux, " " => "" );
                end
                D[ i[ 1 ] ] = aux;
            end
        end
    end
    delete!( Variables, "ExperimentSettings" );
    Variables = merge( Variables, D );
    return Variables
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
    ExperimentDate2String( Variables::Dict{Any, Any} ) -> Variables::Dict{Any, Any}
        Extracting the date of the BRW creation
        using Dates
"""
function ExperimentDate2String( Variables::Dict{Any, Any} )
    X = Variables[ "ExperimentDateTime" ];
    Dt = split( X, ":" );
    Dt[ end ] = string( round( Int, parse( Float64, replace( Dt[ end ], r"[A-Z]" => "" ) ) ) );
    newDt = String( "" );
    for i in Dt
        newDt = string( newDt, ":", i );
    end
    newDt = newDt[ 2:end ];
    X = Dates.DateTime( newDt );
    Variables[ "ExperimentDateTime" ] = string( X );
    T = Dates.format( X, RFC1123Format );
    println( "Creation Date: ", T )
    return Variables
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
    GetSizeChunk(  Variables::Dict{ Any, Any }, limLow::Float64, limHigh::Float64 ) -> ε, ω
        using AccuraTools.size_segments, AccuraTools.div_ab
"""
function GetSizeChunk(  Variables::Dict{ Any, Any }, limLow::Float64, limHigh::Float64 )
    NRecFrames = Variables[ "NRecFrames" ];
    temporal = Variables[ "BRWsize" ] ./ div_ab( NRecFrames );
    SamplingRate = Variables[ "SamplingRate" ];
    aux = findall( temporal .> limLow .&& temporal .< limHigh );
    if isempty( aux )
        ε, ω = size_segments( NRecFrames, SamplingRate );
    else
        n = div_ab( NRecFrames )[ aux[ 1 ] ];
        ω = ceil( Int, ( NRecFrames / n ) ); # number of final frames (chunk size in frames)
        ε = Array{ Int64 }( undef, n, 2 ); # preallocate
        ε[ :, 1 ] = collect( 1:ω:NRecFrames ); # start and
        ε[ :, 2 ] = ε[ :, 1 ] .+ ω .- 1; # end in frames of each chunk (to cut)
        if !isinteger( NRecFrames / n )
            ε[ end, 2 ] = NRecFrames; # end in frames of each chunk (to cut)
        end
    end
    return ε, ω
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
     Get_Chunks( Variables::Dict{ Any, Any }, Output_Chunks::String, limLow = 0.1, limHigh = 0.4 )
        -> Variables::Dict{ Any, Any }
        Cuts the brw dataset into several more manageable segments
        using AccuraTools.GetSizeChunk, JLD, HDF5
"""
function Get_Chunks( Variables::Dict{ Any, Any }, Output_Chunks::String, limLow = 0.1, limHigh = 0.4 )
    Raw = Variables[ "Raw" ];
    Σ = h5read( Variables[ "BRWname" ], Raw );
    TotalRecFrames = size( Σ, 1 );
    nChs = Variables[ "nChs" ];
    NRecFrames = 0;
    try
       NRecFrames = Variables[ "NRecFrames" ];
    catch e
        try;
            NRecFrames = Int( TotalRecFrames/nChs );
            catch e
        end
        if !( maximum( Variables[ "TOC" ] ) == NRecFrames )
            nChs = Int( TotalRecFrames / maximum( Variables[ "TOC" ] ) );
            Variables[ "nChs" ] = nChs;
        end
    end
    Variables[ "NRecFrames" ] = NRecFrames;
    SamplingRate = Variables[ "SamplingRate" ];
    ε, ω = GetSizeChunk(  Variables, limLow, limHigh )
    # number of spaces that occupies the number of seconds of duration of the experiment
    n_char = length( string( Int( floor( NRecFrames / SamplingRate ) ) ) );
    # number of bins
    n = size( ε, 1 );
    # number of spaces occupied by the number of chunks
    n_bins = length( string( n ) );
    #
    for i = 1:n
        # number of B to cut ( 1->4096, 1->ω )
        BIN = Array{ UInt16 }( undef, nChs, ω ); # preallocate
        # values corresponding to the specific BIN. Channel 1,1 has the frames 1, 4097, 8193...etc
        β = collect( ( ε[ i, 1 ] - 1 ):ε[ i, 2 ] );
        finale = ω
        try
             for j = 1:ω
                 # take those frames out of the vector Σ, and put them in array in BIN
                 BIN[ :, j ] = Σ[ ( β[ j ] * nChs )+ 1:( nChs * β[ j + 1 ] ) ];
                 finale = j
             end
        catch e
        end
        if finale != ω
            BIN = BIN[ :, 1:finale ];
        end
        ini = lpad( string( Int(floor( ( β[ 1 ] )/SamplingRate ) ), "s" ), n_char + 1, "0" );
        eni = lpad( string( Int( floor( ( β[ end ] + 1 ) / SamplingRate ) ), "s" ), n_char + 1, "0" );
        bin_time = string( ini, "-", eni );
        BINname = string( "BIN", lpad( i, n_bins, "0" ), "_", bin_time, ".jld" );
        BINname = joinpath( Output_Chunks, BINname );
        save( BINname, "data", BIN );
    end
    return Variables
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
    size_segments( NRecFrames::Int64, SamplingRate::Float64 ) -> ε::Vector, ω::Inf64
        Determines the ideal size of each chunck
        using AccuraTools.div_ab
"""
function size_segments( NRecFrames::Int64, SamplingRate::Float64 )
    if isinteger( SamplingRate ) # thus, there are chunks of 1 second
        if isinteger( NRecFrames/SamplingRate )
            n = Int.( NRecFrames/SamplingRate );
        end
    elseif isinteger( NRecFrames/floor( SamplingRate ) )
        n = Int.( NRecFrames/floor( SamplingRate ) );
    elseif isinteger( NRecFrames/ceil( SamplingRate ) )
        n = Int.( NRecFrames/ceil( SamplingRate ) );
    else
        div_T = div_ab( NRecFrames );
        div_sec = div_T/SamplingRate;
        # range of number of chunks to be cut, normally high to work at ease
        hi = 4; lo = 2; # segundos
        if !isempty( div_T )
            # finds one of the n-frames dividers within the range
            selected_divs = div_T[ findall( hi .>= div_sec .>= lo ) ];
        end
        if !isempty( selected_divs ) # if there was, grab the first one
            n = Int( NRecFrames/selected_divs[ 1 ] );
        else # if there wasn't, a default one.
            n = 60;
        end
    end
    ω = ceil( Int, ( NRecFrames / n ) ); # number of final frames (chunk size in frames)
    ε = Array{ Int64 }( undef, n, 2 ); # preallocate
    ε[ :, 1 ] = collect( 1:ω:NRecFrames ); # start and
    ε[ :, 2 ] = ε[ :, 1 ] .+ ω .- 1; # end in frames of each chunk (to cut)
    if !isinteger( NRecFrames / n )
        ε[ end, 2 ] = NRecFrames; # end in frames of each chunk (to cut)
    end
    return ε, ω
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
    div_ab( n::Int64, lo::Int = 1, hi::Int = n ) -> σ::Vector{Int64}
        Divisors of the number n between the values "lo" and "hi".
        If they are not defined then it takes from 1 to n
"""
function div_ab( n::Int, lo::Int = 1, hi::Int = n )
    ρ = collect( 1:floor( Int, sqrt( n ) ) ) ; # the numbers one by one, from the square root
    σ1 = findall( n.%ρ .== 0 ); # square root divisors ( remainder = 0 )
    σ2 = Int.( ( n ) ./ ( σ1 ) ); # Take out the pairs (of 100, 2-50, 10-10, etc.)
    σ = sort( unique( vcat( σ1, σ2 ) ) ); # remove duplicates, concatenate, sort
    aux1 = @isdefined lo;
    aux2 = @isdefined hi;
    if aux1 && aux2
        rn = σ[ findall( hi .>= σ .>= lo ) ];
        if isempty( rn )
            println(" There is no divisors of $n between $lo and $hi" )
        else
            return rn
        end
    else
        return σ
    end
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #

""" 
convert: Chloe's attempt to update this for the latest documented 3 brain systems.... Chloe voltage conversion new documentation from Chiara @3Brain
"""
function VoltageConvert( data::Matrix{UInt16}, Variables::Dict{Any, Any} )
    MaxVolt = Variables[ "MaxAnalogValue" ];
    MinVolt = Variables[ "MinAnalogValue" ];
    MaxDig = Variables[ "MaxDigitalValue" ];
    MinDig = Variables[ "MinDigitalValue" ];
    BIN = @. MinVolt + (data * (MaxVolt - MinVolt)/ (MaxDig - MinDig));
    return BIN
end
#•·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
    MeanΔxCI( W::Vector, percent::Float64 ) -> xmean, Δx, C1, C2 (::Float64)
    Assumes a Normal distribution. Obtains the confidence interval with "percent" quantile,
    Returns mean, standard error, CI superior and CI inferior
    using BinningAnalysis, Distributions
"""
function MeanΔxCI( W::Vector, percent::Float64 )
    xmean, Δx = BinningAnalysis.jackknife( identity, W );
    d = Distributions.Normal( xmean, std( W ) );
    C1 = xmean + ( Δx * Distributions.quantile( d, percent ) );
    C2 = xmean - ( Δx * Distributions.quantile( d, percent ) );
    return xmean, Δx, C1, C2
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
    StaticThresholdEvents( channel::Vector{Float64}, thr::Real, parameters::Dict{String, Int64} )
    -> index::Vector{Int64}
    Index of detected suprathreshold events with a pre-established threshold
    using AccuraTools.STExbin
"""
function StaticThresholdEvents( channel::Vector{Float64}, thr::Real, parameters::Dict{String, Int64} )
    index = STExbin( channel, thr, parameters );
    return index
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
    Thresholding( W::VecOrMat ) -> Todos::Matrix{Int64}, T::Vector{Float64}
    using HistogramThresholding, ImageContrastAdjustment, Suppressor
"""
function Thresholding( W::VecOrMat )
    W = Float64.( vec( W ) );
    n = length( W );
    edges, conteo = HistogramThresholding.build_histogram( W, length( keys( countmap( W ) ) ) );
    t = zeros( 9 ); Todos = zeros( n, length( t ) );
    @suppress begin
        thr = find_threshold( UnimodalRosin( ), conteo[ 1:end ], edges );
        Todos[ W .>= thr, 1 ] .= 1; t[ 1 ] = thr;
        thr = find_threshold( MinimumIntermodes( ), conteo[ 1:end ], edges );
        Todos[ W .>= thr, 2 ] .= 1; t[ 2 ] = thr;
        thr = find_threshold( Intermodes( ), conteo[ 1:end ], edges );
        Todos[ W .>= thr, 3 ] .= 1; t[ 3 ] = thr;
        thr = find_threshold( MinimumError( ), conteo[ 1:end ], edges );
        Todos[ W .>= thr, 4 ] .= 1; t[ 4 ] = thr;
        thr = find_threshold( Moments( ), conteo[ 1:end ], edges );
        Todos[ W .>= thr, 5 ] .= 1; t[ 5 ] = thr;
        thr = find_threshold( Otsu( ), conteo[ 1:end ], edges );
        Todos[ W .>= thr, 6 ] .= 1; t[ 6 ] = thr;
        thr = find_threshold( Entropy( ), conteo[ 1:end ], edges );
        Todos[ W .>= thr, 7 ] .= 1; t[ 7 ] = thr;
        thr = find_threshold( Balanced( ), conteo[ 1:end ], edges );
        Todos[ W .>= thr, 8 ] .= 1; t[ 8 ] = thr;
        thr = find_threshold( Yen( ), conteo[ 1:end ], edges );
        Todos[ W .>= thr, 9 ] .= 1; t[ 9 ] = thr;
    end
    return Int.( Todos ), t
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
    ThresholdingPlots( Todos::Matrix{Int64} ) -> TF::Plot
    using AccuraTools.Zplot, Plots
"""
function ThresholdingPlots( Todos::Matrix{Float64} )
    xmessage = string( "UnimodalRosin : ", length( findall( Todos[ :, 1 ] .== 1 ) )," channels" );
    T1 = Zplot( Todos[ :, 1 ], "W" ); T1 = plot!( title = xmessage, cbar = :none );
    xmessage = string( "MinimumIntermodes: ", length( findall( Todos[ :, 2 ] .== 1 ) ),"channels" );
    T2 = Zplot( Todos[ :, 2 ], "W" ); T2 = plot!( title = xmessage, cbar = :none );
    xmessage = string( "Intermodes: ", length( findall( Todos[ :, 3 ] .== 1 ) )," channels" );
    T3 = Zplot( Todos[ :, 3 ], "W" ); T3 = plot!( title = xmessage, cbar = :none );
    xmessage = string( "MinimumError: ", length( findall( Todos[ :, 4 ] .== 1 ) )," channels" );
    T4 = Zplot( Todos[ :, 4 ], "W" ); T4 = plot!( title = xmessage, cbar = :none );
    xmessage = string( "Moments: ", length( findall( Todos[ :, 5 ] .== 1 ) )," channels" );
    T5 = Zplot( Todos[ :, 5 ], "W" ); T5 = plot!( title = xmessage, cbar = :none );
    xmessage = string( "Otsu: ", length( findall( Todos[ :, 6 ] .== 1 ) )," channels" );
    T6 = Zplot( Todos[ :, 6 ], "W" ); T6 = plot!( title = xmessage, cbar = :none );
    xmessage = string( "Entropy: ", length( findall( Todos[ :, 7 ] .== 1 ) )," channels" );
    T7 = Zplot( Todos[ :, 7 ], "W" ); T7 = plot!( title = xmessage, cbar = :none );
    xmessage = string( "Balanced: ", length( findall( Todos[ :, 8 ] .== 1 ) )," channels" );
    T8 = Zplot( Todos[ :, 8 ], "W" ); T8 = plot!( title = xmessage, cbar = :none );
    xmessage = string( "Yen: ", length( findall( Todos[ :, 9 ] .== 1 ) )," channels" );
    T9 = Zplot( Todos[ :, 9 ], "W" ); T9 = plot!( title = xmessage, cbar = :none );

    TF = plot(
        T1, T2, T3, T4, T5, T6, T7, T8, T9,
        layout = ( 3, 3 ), wsize = ( 700, 700 ), titlefont = ( 8, "arial" ) );
    return TF
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
    Z0( X::VecOrMat, nChs::Int64 ) -> Z::Matrix{Int64}
    using Plots
"""
function Z0( X::VecOrMat, nChs::Int64 )
    X = vec( X );
    Z = zeros( Int, nChs );
    n = Int( sqrt( nChs ) );
    Z[ X ] = Z[ X ] .+ 1;
    Z = reverse( reshape( Z, n, n )', dims = 1 );
    return Z
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
    ZW( X::VecOrMat ) -> Z::Matrix{typeof(X)}
    using Plots
"""
function ZW( X::VecOrMat )
    X = vec( X );
    n = Int( sqrt( length( X ) ) );
    Z = reverse( reshape( X, n, n )', dims = 1 );
    return Z
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
    Zplot( Z::Matrix, which::String, cm = :greys, nChs = 4096 ) -> F::Plot
    using Plots, AccuraTools.Z0, AccuraTools.ZW
"""
function Zplot( Z::VecOrMat, which::String, cm = :greys, nChs = 4096 )
    if which == "0"; Z = Z0( Z, nChs ); elseif which == "W"; Z = ZW( Z ); end
    F = heatmap( Z, aspect_ratio = 1, c = cm, axis = ( [ ], false ), wsize = ( 400, 400 ) );
    return F
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #

  #Function to detect positive thresholds... This is the equivalent of the STExbin function by Isabel. 
"""
    PosT( bin::Vector{Float64}, posThr::Real, parameters::Dict{String, Int64} ) -> Index::Vector{Int64}
"""
function PosT( bin::Vector{Float64}, posThr::Real, parameters::Dict{String, Int64} ) #is this correct? 
    
   #posThr = abs(thr);
    distance = parameters["distance"];
    PTX =findall( bin .>= posThr); # positive threshold crossings (X)
    part_index = [ ];
    if !isempty( PTX )
       b = 1;
        while b == 1
             mates = diff( PTX ); # distancia entre ellos
            nachbar = findall( mates .<= distance ) .+ 1; # cuales estan cerca
            if isempty( nachbar )
                b = 0;
            else 
                removethem = zeros( Int, size( nachbar, 1 ));
    for i = 1:size(nachbar, 1 ) #this bit looks for where the peak is
                    if !isless((bin[ PTX[ nachbar[ i ] ] ]), ( bin[ PTX[ nachbar[ i ] - 1 ]] ) ) 
                        removethem[ i ] = nachbar[ i ] - 1;
                    else #i.e. now we're at the peak
                        removethem[ i ] = nachbar[ i ];
                    end
                end
                if size( removethem, 1 ) > 1
                    PTX[ unique(removethem)].= 0;
                    filter!(x -> x !=0, PTX);
                    else 
                    PTX = PTX[ Bool.( PTX .!= PTX[ removethem[ 1 ] ] ) ];
                end
            end
        end
        push!( part_index, PTX)
    else
        push!( part_index, [] )
    end
    return sort( unique( vcat(part_index...) ) )
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
#now to find the pos events.... 
"""
  PTXchannel( channel::Vector{Float64}, parameters::Dict{String, Int64} )
    -> thrs::Vector{Float64}, real_index::Vector{Inf64}
    using AccuraTools.donoho, AccuraTools.PosT
"""
function PTXchannel( channel::Vector{Float64}, parameters::Dict{String, Int64} )
    window = parameters[ "window" ];
    bit = parameters[ "bit" ];
    σ = parameters[ "cte" ];
    real_index = [ ];
    i = 1;
    M = Int( ( ( ( i - 1 ) * bit ) + 1 ) ); N = Int( M + window - 1 );
    posThr = [ ];
    thr2 = 0;
    while N <= ( length( channel ) - Int( window - 1 ) )
        bin = channel[ M:N ];
        if !iszero( bin ) 
            thr2 = σ*( donoho( bin ) );
            part_index = PosT( bin, thr2, parameters );
            if !isempty( part_index)
                real_index = vcat( real_index, ( part_index .+ M .- 1 ) );
            end
        end
        i = i + 1;
        M = Int( ( ( ( i - 1 ) * bit ) + 1 ) ); N = Int( M + window - 1 );
    end
    real_index = unique( real_index );
    if !isempty( real_index )
        push!( posThr, thr2 )
    else
        push!( posThr, [0] )
    end
    return vcat( posThr... ), real_index
end
# •·•·•·•·•·
"""
    STExbin( bin::Vector{Float64}, thr::Real, parameters::Dict{String, Int64} ) -> Index::Vector{Int64}
"""
function STExbin( bin::Vector{Float64}, thr::Real, parameters::Dict{String, Int64} )
    distance = parameters[ "distance" ];
    ST = findall( bin .<= thr ); # eventos que pasan el umbral
    index_parcial = [ ];
    if !isempty( ST )
        a = 1;
        while a == 1
            distances = diff( ST ); # distancia entre ellos
            nears = findall( distances .<= distance ) .+ 1; # cuales estan cerca
            if isempty( nears )
                a = 0;
            else
                remove = zeros( Int, size( nears, 1 ) );
                for i = 1:size( nears, 1 )
                    # si el primero es menor que el segundo, quita el segundo
                    if isless( ( bin[ ST[ nears[ i ] ] ] ),
                            ( bin[ ST[ nears[ i ] - 1 ] ] ) )
                        remove[ i ] = nears[ i ] - 1;
                    else
                        remove[ i ] = nears[ i ];
                    end
                end
                if size( remove, 1 ) > 1
                    ST[ unique( remove ) ] .= 0;
                    filter!( x -> x != 0, ST );
                else
                    ST = ST[ Bool.( ST .!= ST[ remove[ 1 ] ] ) ];
                end
            end
        end
        push!( index_parcial, ST )
    else
        push!( index_parcial, [ ] )
    end
    return sort( unique( vcat( index_parcial... ) ) )
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
    STExChannel( channel::Vector{Float64}, parameters::Dict{String, Int64} )
    -> thrs::Vector{Float64}, index_real::Vector{Inf64}
    using AccuraTools.donoho, AccuraTools.STExbin
"""
function STExChannel( channel::Vector{Float64}, parameters::Dict{String, Int64} )
    window = parameters[ "window" ];
    bit = parameters[ "bit" ];
    σ = parameters[ "cte" ];
    index_real = [ ];
    i = 1;
    I = Int( ( ( ( i - 1 ) * bit ) + 1 ) ); J = Int( I + window - 1 );
    thrs = [ ];
    thr = 0;
    posthr = 0;
    while J <= ( length( channel ) - Int( window - 1 ) )
        bin = channel[ I:J ];
        if !iszero( bin ) 
            thr = -1*σ*abs( donoho( bin ) );
            index_parcial = STExbin( bin, thr, parameters );
            if !isempty( index_parcial )
                index_real = vcat( index_real, ( index_parcial .+ I .- 1 ) );
            end
        end
        i = i + 1;
        I = Int( ( ( ( i - 1 ) * bit ) + 1 ) ); J = Int( I + window - 1 );
    end
    index_real = unique( index_real );
    if !isempty( index_real )
        push!( thrs, thr )
    else
        push!( thrs, [0] )
    end
    return vcat( thrs... ), index_real
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
    Noise-adaptive Optimal Thresholding
        donoho( x::Vector ) -> thr::Float64
        using Statistics
https://www.nature.com/articles/s41598-021-93088-w#Sec2
"""
@inline donoho( x ) = ( median( abs.( x ) ) / 0.6745 );
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
    FiltMUAremez( Variables::Dict{Any, Any}, channel::Vector{Float64} -> MUA::Vector{Float64}
    using DSP
"""
function FiltMUAremez( Variables::Dict{Any, Any}, cortex_ch::Vector{Int64} )
    lF = 300;
    fac = 10;
    HF = 3000;
    SamplingRate = Variables[ "SamplingRate" ];
    NYQ = floor( Int, SamplingRate / 2 );
    order = Int( floor( ( SamplingRate / lF ) / 5 ) );
    bpass = remez(
        ( order + 1 ), [ ( 0, lF - fac ) => 0, ( lF, HF ) => 1, ( HF + fac, NYQ ) => 0 ],
            Hz = SamplingRate );
    MUA = filtfilt( DSP.Filters.PolynomialRatio( bpass, [ 1.0 ] ), cortex_ch );
    return MUA
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
    FilteringMatrix( data::Matrix{Float64}, Variables::Dict{Any, Any} ) -> datafilt::Matrix{Float64}
    using AccuraTools.FiltroMUAremez
"""
function FilteringMatrix( datactx::Matrix{Float64}, Variables::Dict{Any, Any} )
    datafilt = copy( datactx );
    for k = 1:size( datactx, 1 )
        channel = datactx[ k, : ];
        MUA = FiltMUAremez( Variables, channel );
        datafilt[ k, : ] = MUA;
    end
    return datafilt
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"Try to convert the ^^ function to Dict String, Any
"
function FiltM( datactx::Matrix{Float64}, Variables::Dict{Any, Any})
    datafilt = copy( datactx );
    for k = 1:size( datactx, 1 )
        channel = datactx[ k, : ];
        MUA = FiltMUAremez( Variables, cortex_ch );
        datafilt[ k, : ] = MUA;
    end
    return datafilt
end
"""
    ms2frames( time::Real, Variables::Dict{Any, Any} ) -> x::Float64
"""
function ms2frames( time::Real, Variables::Dict{Any, Any} )
    SamplingRate = Variables[ "SamplingRate" ];
    if time != 0; x = ceil( Int, ( time * SamplingRate ) / 1000 ); else; x = 1; end
    return x
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
    ΔV( Variables::Dict{Any, Any}, BIN::Matrix{Float64}, ΔT::Int64 ) -> STD::Vector{Float64}
        using AccuraTools.ms2frames
"""
function ΔV( Variables::Dict{Any, Any}, BIN::Matrix{Float64}, ΔT::Int64 )
    ΔT = ms2frames( ΔT, Variables );
    STD = vec( std( ( BIN - circshift( BIN, ( 0, ΔT )  ) ), dims = 2 ) );
    return STD
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
    searchdir( path::String, key::String ) -> Vector{String}
    Find the files with the key word on their name in the given path, and returns a vector of strings with the
    fullname of those files (complete path)
"""
@inline searchdir( path::String, K::String ) = filter( x -> endswith( x, K ), readdir( path; join = true ) );
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
    PositiveSaturation( data::Matrix{UInt16}, Variables::Dict{Any, Any} )
        -> DiscartedChannels::Vector{Int64}, earth::Int64, data::Matrix{Float64}
    Patching the frames with maximum voltage value permitted with random values from the same channel.
    If the symmetry of the values of each channel is equal to 0, it is a member of the discarded list of
    channels or is the earth. Naturally, is expected that only one channel is the reference channel.
    using StatsBase
"
#function VoltegeConversion( data::Matrix{UInt16}, Variables::Dict{Any, Any} )
#    MaxVolt = Variables[ "MaxAnalogValue" ];
#    MinVolt = Variables[ "MinAnalogValue" ];
 #   SignalInversion = Variables[ "ScaleFactor" ];
 #   BitDepth = 12;
  #  ADCCountsToMV = SignalInversion * ( ( MaxVolt - MinVolt )/ ( 2 ^ BitDepth ) );
 #   MVOffset = SignalInversion * MinVolt;
 #   BIN = @. MVOffset + ( data * ADCCountsToMV );
#    return BIN
#end"""

"""#from Isabel's github :-> 
function SaturacionPositiva( Variables::Dict{String, Any}, data )

    nChs = length( Variables[ "StoredChIdxs" ] );
    SamplingRate = Variables[ "SamplingRate" ];
    MinVolt = Float64( Variables[ "MinAnalogValue" ] + 1 );
    channels = collect( 1:nChs );
    binsize = size( data, 2 );
    MinVoltSat = sum( Int.( data .== MinVolt ), dims = 2 );

    if length( findall( vec( MinVoltSat ) .== 0 ) ) == ( nChs - 1 )
        earth = setdiff( channels, findall( vec( MinVoltSat ) .== 0 ) )[ ];
        MinVoltSatChannels = [ ];
        data[ earth, : ] .= 0;
    else
        SatFramesLimit = ceil( Int, SamplingRate*0.005 );
        NineNinePercent = ceil( Int, binsize*0.99 );
        MinVoltSatChannels = findall( vec( NineNinePercent .>= MinVoltSat .>= SatFramesLimit ) );
        aux = findall( vec( MinVoltSat ) .>= NineNinePercent );
        if length( aux ) == 1
            earth = aux[ 1 ];
            data[ earth, : ] .= 0;
        else
            println( "there is something fishy here" );
            earth = []
        end
    end
    return data, earth, MinVoltSatChannels
end"""

function PositiveSaturation( data::Matrix{Float64}, Variables::Dict{Any, Any} )
    MaxAV = Float64( Variables[ "MaxAnalogValue" ] - 10 );
    SamplingRate = Variables[ "SamplingRate" ];
    binsize = size( data, 2 );
    MaxVoltSatChannels = findall( vec( sum( Int.( data .>= MaxAV ), dims = 2 ) .!= 0 ) )
    MaxVoltFrames = [ ];
    for k = 1:length( MaxVoltSatChannels )
        push!( MaxVoltFrames, findall( data[ MaxVoltSatChannels[ k ], : ] .>= MaxAV ) );
    end
    for i = 1:length( MaxVoltSatChannels )
        channel = MaxVoltSatChannels[ i ];
        nframes = length( MaxVoltFrames[ i ] );
        valid = data[ channel, setdiff( 1:binsize, MaxVoltFrames[ i ] ) ];
        patch = sample( valid, nframes );
        data[ channel, MaxVoltFrames[ i ] ] = patch;
    end
   # DiscartedChannels = findall( vec( sum( data, dims = 2 ) ) .== 0 );
    return data
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
    Get_Groups( W::Vector{Int64} ) -> grupos::Vector{Any}, loose::( Vector{Any} )
       Groups of adjacent channels are formed from the initial indexes W. Those channels separated from the bulk
        of the others are considered loose ones.
        using AccuraTools.neighborgs, AccuraTools.reverberation
"""
function Get_Groups( W::Vector{Int64} )
    grupos = [ ];
    for i in W
        _, vecinos = neighborgs( i, 1 );
        grupo = sort( intersect( vecinos, W ) );
        if isempty( grupo )
            push!( grupos, i )
        else
            push!( grupos, vcat( grupo, i ) )
        end
    end
    loose = grupos[ findall( length.( grupos ) .== 1 ) ];
    deleteat!( grupos, findall( length.( grupos ) .== 1 ) );
    a = length( grupos ); grupos = reverberation( grupos ); b = length( grupos );
    while a != b
        a = length( grupos ); grupos = reverberation( grupos ); b = length( grupos );
    end
    return grupos, loose
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
    reverberation( grupos::Vector{Any} ) -> more_groups::Vector{Any}
        Intermediate step for Get_Groups function
"""
function reverberation( grupos::Vector{Any} )
    adjoining_channels = 0
    for i in grupos
        adjoining_channels = vcat( adjoining_channels, i )
    end
    adjoining_channels = adjoining_channels[ adjoining_channels .!= 0 ];
    adjoining_channels = unique( adjoining_channels );
    more_groups = [ ]
    for i in adjoining_channels
        temporal = 0
        for j = 1:length( grupos )
            if !isempty( findall( grupos[ j ] .== i ) )
                temporal = vcat( temporal, grupos[ j ] );
            end
        end
        temporal = temporal[ temporal .!= 0 ];
        new_group = sort( unique( temporal ) );

        if !isempty( new_group )
            if !isempty( more_groups )
                temp = last( more_groups );
                if !isequal( temp, new_group )
                    push!( more_groups, new_group );
                end
            else
                push!( more_groups, new_group );
            end
        end
    end
    return more_groups
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
    neighborgs( C::Int64, d::Int64 ) ->
        -> A = Array( ( d*2 ) + 1, ( d * 2 ) + 1 ), v = vec( 2*( ( d * 2 ) + 1 ) - 1 );
        The d-neighborhood is calculated from the channel (C) as a center
        A = array where C is the center and is in chip order
        v = same neighboring channels as A but in vector form and without C ( 8 channels )
"""
function neighborgs( C::Int64, d::Int64 )
    Layout = reverse( reshape( collect( 1:4096 ), 64, 64 )', dims = 1 );
    x_c = findall( Layout .== C )[ ][ 2 ]; y_c = findall( Layout .== C )[ ][ 1 ];
    aux = [ ( x_c - d ),( x_c + d ), ( y_c - d ), ( y_c + d ) ]
    aux[ aux .< 1 ] .= 1; aux[ aux .> 64 ] .= 64;
    A = Layout[ aux[ 3 ]:aux[ 4 ], aux[ 1 ]:aux[ 2 ] ];
    v = vec( A )[ vec( A ) .!= C ];
    return A, v
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
function FigureGroups( grupos::Vector, loose::Vector = [ ], cm = :twilight )
    Z = zeros( Int, 4096 );
    for i = 1:size( grupos, 1 )
        Z[ grupos[ i ] ] .= floor( Int, log( length( grupos[ i ] ) ) ) + 2 ;
    end
    Z[ Int.( loose ) ] .= 1
    Z = reverse( reshape( Z, 64, 64 )', dims = 1 )
    F = heatmap(
        Z,
        aspect_ratio = 1,
        c = cm,
        axis = ( [ ], false ),
        wsize = ( 400, 400 ),
        cbar = :none );
    return F
end

# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
    FillingHolesCrux( Seleccion::Vector{Int64} ) -> NewSelection::Vector{Int64}
        using AccuraTools.neighborgs
"""
function FillingHolesCrux( Seleccion::Vector{Int64} )

    Z = reverse( reshape( 1:4096, 64, 64 )', dims = 1 );
    arriba = Z[ 64, 2:( end - 1 ) ]; abajo = Z[ 1, 2:( end - 1 ) ];
    izquierda = Z[ 2:( end - 1 ), 1 ]; derecha = Z[ 2:( end - 1 ), 64 ];
    esquinas = vcat( Z[ 1, 1 ], Z[ 64, 64 ], Z[ 64, 1 ], Z[ 1, 64 ] );
    bordes = vcat( arriba, abajo, izquierda, derecha, esquinas );
    bordes_noesquinas = vcat( arriba, abajo, izquierda, derecha );

    P = true
    while P
        X1 = setdiff( setdiff( setdiff( 1:4096, Seleccion ) ), bordes );
        p1 = [ ];
        for i = 1:length( X1 )

            x = X1[ i ];
            A, _ = neighborgs( x, 1 );
            CruzH = [ A[ 2, 1 ], A[ 2, 3 ] ];
            CruzV = [ A[ 1, 2 ], A[ 3, 2 ] ];

            x1 = length( CruzH[ CruzH .∈ [ Seleccion ] ] ) == 2;
            x2 = length( CruzV[ CruzV .∈ [ Seleccion ] ] ) == 2;

            if ( x1 || x2 )
                push!( p1, x )
            end
        end

        X1 = intersect( arriba, setdiff( setdiff( 1:4096, Seleccion ) ) );

        p2 = [ ];
        for i = 1:length( X1 )

            x = X1[ i ];
            A, _ = neighborgs( x, 1 );
                Cruz = [ A[ 2, 1 ], A[ 2, 3 ] ];
                x1 = length( Cruz[ Cruz .∈ [ Seleccion ] ] ) == 2;
            if x1
                push!( p2, x )
            end
        end

        X1 = intersect( abajo, setdiff( setdiff( 1:4096, Seleccion ) ) );

        p3 = [ ];
        for i = 1:length( X1 )

            x = X1[ i ];
            A, _ = neighborgs( x, 1 );
                Cruz = [ A[ 1, 1 ], A[ 1, 3 ] ];
                x1 = length( Cruz[ Cruz .∈ [ Seleccion ] ] ) == 2;
            if x1
                push!( p3, x )
            end
        end

        X1 = intersect( vcat( izquierda, derecha ), setdiff( setdiff( 1:4096, Seleccion ) ) );

        p4 = [ ];
        for i = 1:length( X1 )

            x = X1[ i ];
            A, _ = neighborgs( x, 1 );
                Cruz = [ A[ 1, 2 ], A[ 3, 2 ] ];
                x1 = length( Cruz[ Cruz .∈ [ Seleccion ] ] ) == 2;
            if x1
                push!( p4, x )
            end
        end

         X1 = intersect( esquinas, setdiff( setdiff( 1:4096, Seleccion ) ) );
         p5 = [ ];
         for i = 1:length( X1 )
             x = X1[ i ];
             _, v = neighborgs( x, 1 );

             if length( v[ v .∈ [ Seleccion ] ] ) >= 2
                 push!( p5, x )
             end
         end
        posibles = [ ];
        posibles = vcat( Int.( p1 ), Int.( p2 ), Int.( p3 ), Int.( p4 ), Int.( p5 ) );

        if isempty( posibles )
            P = false
        else
            Seleccion = vcat( posibles, Seleccion );
        end

    end
    return Int.( Seleccion )
end
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
    Cut_Spikes( index_final::Vector{Any}, MUA::Vector{Float64}, parametros::Dict{String, Int64} )
        => spikes::Matrix{Float64}
        takes the centers detected in the ChannelXEvents step and cuts the spikes, with w_pre ms before
          and w_pos later.

    """
function Cut_Spikes( index_final::Vector{Int64}, MUA::Vector{Float64}, parametros::Dict{String, Int64} )
    w_post = parametros[ "w_post" ]; w_pre = parametros[ "w_pre" ];
    index_final = index_final[ Bool.( 1 .- ( ( index_final .+ w_post ) .> length( MUA ) ) ) ];
    index_final = index_final[ Bool.( 1 .- ( ( index_final .- w_pre ) .< 1 ) ) ];
    spikes = zeros( length( index_final ), ( w_pre + w_post + 1) );
    for i = 1:length( index_final )
        spikes[ i, : ] = MUA[ ( index_final .- w_pre )[ i ]:( index_final .+ w_post )[ i ] ];
    end
    return spikes
end

end  #•·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·• #
