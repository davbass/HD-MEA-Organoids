"""
Tools to find and delete recalibration artefacts, from 3Brain STIMULO Chips, with ~11000 Hz sampling rate. Works with BRW (HDF5) files, from Brainwave 5 program. 
Work developed by Dr. Chloe Hall, from AG Mittmann at UM Mainz. 
hallchlo@uni-mainz.de
https://github.com/chloemhall
""" 
module Artefacts

using HDF5
using JLD 
using Plots
using Statistics
using StatsBase
using DelimitedFiles
using FFTW
using DSP
using JSON
#using SignalProcessing

export findTimestamps;
export realNeighbours;
export extract_channel_indexes; 
export apply_highpass_filter;
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
Find the timestamps of recalibration, which is where we need to delete. Adds an extra 20ms gap, to create a barrier i.e. 200 pts- 100 pts either side of the recalibration timestamp. Returns a matrix named "todelete " which will cause the 
"""
function findTimestamps( data::Matrix{Float64} ) # not sure if this is correct.. 

calCh = data[1, :]; # input here, which channel is the Earth 
cutoff = 3;
indexArray = findall(x -> x > cutoff, calCh); # finds all the points 
todelete = [ ];
    mat2 = [ ];
    mat3 = [ ];
barrier = collect(1:100); #creates a barrier of points to delete.
    # now we're trying to make sure we only delete the correct points...
    for i = 3:length(indexArray)
        
        if indexArray[i]-indexArray[i-2] == 2
        addpos = (indexArray[i] .+ barrier);
        addneg = (indexArray[i-2] .- barrier)
        mat2 = push!(mat2, addpos );
        mat3 = push!(mat2, addneg );
        else continue
        end
        
    end
    

#^^  these two are the correct data points, that I want to now delete from data 
mat22 = collect(Iterators.flatten(mat2)); #this flattens the matrices I created, into one continuous vector.
mat33 = collect(Iterators.flatten(mat3));
todelete = [mat22; mat33];
return todelete

end 

# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
"""
 function to make a nested vector of all channels in SS cortex and their neighbours which were also selected. 
"""
function realNeighbours( cortex_ch::Vector{Int64} , checkCh::Vector{Any} ) 

neighbourGrid = [ ];
    list=[ ];
for i in 1:length(checkCh)
    list = checkCh[i];
    newB = intersect(cortex_ch, list)
    neighbourGrid = push!(neighbourGrid, newB)
        
end # output is neighbourGrid nested vector. 
    return neighbourGrid
end # ie. function end
#·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
#·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
function extract_channel_indexes(json_file::AbstractString, group_name::AbstractString)
    Chdata = JSON.parsefile(json_file) # data file with the channel selected. 
    groups = Chdata["Groups"]
    
    for group in groups
        if group["UserDefinedName"] == group_name
            channel_indexes = group["PixelIndexes"]
            return channel_indexes
        end
    end
    
    return nothing
end
#·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
function apply_highpass_filter(data::Matrix, sample_rate::Float64, cutoff_freq::Float64)
    num_rows, num_cols = size(data)
    filtered_data = similar(data)
    
    for i in 1:num_rows
        row_data = data[i, :]
        b, a = butter(1, cutoff_freq / (sample_rate / 2), "high")
        filtered_row_data = filtfilt(b, a, row_data)
        filtered_data[i, :] = filtered_row_data
    end
    
    return filtered_data
end

#·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #

end #i.e module end

#·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #