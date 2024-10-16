# SPDX-License-Identifier: BSD-2-Clause

module TraitExtraction

using XLSX, CSV
using DataFrames, Unitful, Statistics, Phylo

function uparse_missing(::Missing)
    return missing
end

function uparse_missing(val::AbstractString)
    return uparse(replace(val, r" " => ""))
end

function parse_missing_range(_, ::Missing)
    return missing
end

function parse_missing_range(T::Type{<:Number}, val::AbstractString)
    vec = split(val, r" *- *")
    if length(vec) == 1
        if occursin(r",", val)
            return parse_missing(Vector{T}, vec[1])
        else
            return parse.(T, vec)
        end
    elseif length(vec) == 2
        vv = parse.(T, vec)
        collect(vv[1]:vv[2])
    else
        @warn "Could not parse $val"
        return missing
    end
end

function parse_missing_2n(::Missing)
    return missing
end

function parse_missing_2n(val::AbstractString)
    v2 = replace(val, r"[,.] ?$" => "", r" ?2n ?= ?" => "", r" ?\| ?" => ", ",
                 r" ?, ?" => ", ")
    if !occursin(r"[^ ]", v2)
        return missing
    elseif occursin(r"[^0-9, ]", v2)
        @warn "Could not parse $val"
        return missing
    else
        v = split(v2, r" *, *", keepempty = true)
        vv = tryparse.(Int, v[occursin.(r"[0-9]", v)])
        if any(isnothing, vv)
            @warn "Could not parse $val"
            return missing
        else
            return unique(sort(vv))
        end
    end
end

function parse_missing(_, ::Missing)
    return missing
end

function parse_missing(::Type{Vector{String}}, val::AbstractString)
    v = split(val, r" *, *", keepempty = true)
    v = sort(v[occursin.(r"[0-9a-zA-Z]", v)])
    return isempty(v) ? missing : v
end

function parse_missing(T::Type{<:Vector{<:Number}}, val::AbstractString)
    return sort(parse.(eltype(T), split(val, r" *, *")))
end

function parse_missing(T::Type{<:Number}, val::AbstractString)
    return parse(T, val)
end

function removelabel!(col::Vector{Missing}, _, _)
    return col
end

function removelabel!(col, old::AbstractVector)
    if nonmissingtype(eltype(col)) <: AbstractVector
        for i in eachindex(col)
            if !ismissing(col[i]) && col[i] == old
                col[i] = missing
            end
        end
    else
        @warn "Cannot memove a vector $old from a non-vector column"
    end

    return col
end

function removelabel!(col, old)
    if nonmissingtype(eltype(col)) <: AbstractVector
        for i in eachindex(col)
            ismissing(col[i]) && continue
            if old ∈ col[i]
                filter!(≠(old), col[i])
            end
            if isempty(col[i])
                col[i] = missing
            end
        end
    else
        col[.!ismissing.(col) .& (col .== old)] .= missing
    end

    return col
end

function relabel!(col::Vector{Missing}, ::AbstractString, _)
    return col
end

function relabel!(col, old::AbstractVector, new::AbstractVector)
    if nonmissingtype(eltype(col)) <: AbstractVector
        for i in eachindex(col)
            if !ismissing(col[i]) && old == col[i]
                col[i] = new
            end
        end
    else
        @warn "Cannot relabel a non-vector column from $old to $new"
    end

    return col
end

function relabel!(col, old, new)
    if nonmissingtype(eltype(col)) <: AbstractVector
        for i in eachindex(col)
            if !ismissing(col[i]) && old ∈ col[i]
                col[i] = sort(unique(replace(col[i], old => new)))
            end
        end
    else
        col[.!ismissing.(col) .& (col .== old)] .= new
    end

    return col
end

function clean_traits(df_raw)
    df = copy(df_raw)

    col = df[!, "flower architecture"]
    removelabel!(col, "catkin")
    removelabel!(col, "pentamerous")

    col = df[!, "fruit structure"]
    removelabel!(col, ["dry", "fleshy"])

    col = df[!, "petal fusion"]
    relabel!(col, "fusion", "fused")
    removelabel!(col, ["free", "fused"])

    col = df[!, "life form"]
    removelabel!(col, "hydrophyte")
    removelabel!(col, "geophyte")
    removelabel!(col, ["annual", "perennial"])
    removelabel!(col, ["annual", "biennials", "perennial"])

    col = df[!, "spinescence"]
    removelabel!(col, ["armed", "spinescent", "unarmed"])
    removelabel!(col, ["armed", "spinescent", "spinulate", "unarmed"])
    relabel!(col, ["armed", "spinescent", "spinulate"],
             ["spinescent", "spinulate"])
    relabel!(col, ["armed", "spinescent"], ["spinescent"])
    relabel!(col, ["armed", "spinulate"], ["spinulate"])
    removelabel!(col, ["armed", "unarmed"])
    relabel!(col, ["spinescent", "unarmed"], ["spinescent"])
    relabel!(col, ["spinulate", "unarmed"], ["spinulate"])

    col = df[!, "fruit dehiscence"]
    relabel!(col, ["dehiscent", "indehiscent", "valvate"], ["valvate"])
    relabel!(col, ["dehiscent", "suture"], ["suture"])

    col = df[!, "habitat"]
    relabel!(col, ["aquatic", "aquatic/hygrophilous"], ["hygrophilous"])
    relabel!(col, ["epiphyte", "terrestrial"], ["epiphyte"])

    col = df[!, "succulence"]
    relabel!(col, "flexhy", "fleshy")
    relabel!(col, "yes", "succulent")
    removelabel!(col, ["fleshy", "succulent"])

    col = df[!, "flower sex"]
    relabel!(col, ["bisexual", "pistillate", "staminate"], ["unisexual"])
    relabel!(col, ["bisexual", "pistillate", "unisexual"],
             ["bissexual", "pistillate"])
    relabel!(col, ["pistillate", "staminate", "unisexual"], ["unisexual"])
    relabel!(col, ["staminate", "unisexual"], ["staminate"])
    relabel!(col, ["pistillate", "unisexual"], ["unisexual"])

    col = df[!, "leaf arrangement"]
    relabel!(col, ["bipinnate", "lobed", "pinnate"], ["bipinnate"])
    relabel!(col, ["bipinnate", "pinnate"], ["bipinnate"])
    relabel!(col, ["digitate", "pinnate"], ["digitate"])
    relabel!(col, ["lobed", "pinnate"], ["pinnate"])
    relabel!(col, ["lobed"], ["simple"])
    relabel!(col, ["lobed", "simple"], ["simple"])
    removelabel!(col, ["lobed", "pinnate", "simple"])
    removelabel!(col, ["simple", "trifoliolate"])

    col = df[!, "leaf position"]
    removelabel!(col, ["alternate", "opposite"])
    removelabel!(col, ["alternate", "fascicled", "opposite"])
    relabel!(col, ["decussate", "opposite"], ["decussate"])
    relabel!(col, ["alternate", "spirally"], ["spirally"])
    relabel!(col, "fascicled", "fasciculate")
    relabel!(col, ["verticillate", "whorled"], ["whorled"])
    removelabel!(col, ["alternate", "fasciculate"])

    col = df[!, "fruit shape"]
    relabel!(col, ["cylindrical", "fusiform"], ["fusiform"])
    relabel!(col, ["cylindrical", "ovoid"], ["ovoid"])
    relabel!(col, "turbinate", "ovoid")

    col = df[!, "inflorescence arrangement"]
    relabel!(col, ["clustered", "corymb", "raceme"], ["clustered", "corymb"])
    relabel!(col, ["clustered", "panicle", "raceme"], ["clustered", "panicle"])
    relabel!(col, ["clustered", "panicle", "raceme", "thyrse"],
             ["clustered", "thyrse"])
    relabel!(col, ["raceme", "thyrse"], ["thyrse"])
    relabel!(col, ["raceme", "umbel"], ["umbel"])
    relabel!(col, ["panicle", "thyrse"], ["thyrse"])
    relabel!(col, ["corymb", "panicle", "raceme"], ["corymb", "panicle"])
    relabel!(col, ["panicle", "raceme", "umbel"], ["panicle", "umbel"])

    col = df[!, "habit"]
    for i in 1:nrow(df)
        if ismissing(col[i])
            continue
        end
        if occursin(r"^climber", col[i][1])
            col[i] = ["climber"]
        elseif occursin(r"^character", col[i][1])
            col[i] = ["climber"]
        end
    end
    relabel!(col, ["erect leafy", "herb"], ["erect leafy"])
    relabel!(col, ["erect leafy", "herb", "shrub"], ["erect leafy"])
    relabel!(col, ["herb", "prostrate"], ["prostrate"])
    relabel!(col, ["herb", "prostrate", "shrub"], ["prostrate"])
    relabel!(col, ["herb", "tree", "tree/shrub"], ["herb"])
    relabel!(col, ["tree", "tussock"], ["tree"])
    relabel!(col, ["shrub", "tussock"], ["shrub"])
    relabel!(col, ["prostrate", "shrub", "tree/shrub"], ["prostrate"])
    relabel!(col, ["prostrate", "tussock"], ["prostrate"])
    relabel!(col, ["herb", "tussock"], ["erect leafy"])

    for col in names(df)
        if all(ismissing.(df[!, col]))
            select!(df, Not(col))
            @info "Dropping column $col"
        end
    end

    return df
end

function read_trait_data(filename::AbstractString)
    if occursin(r".xlsx$", filename)
        data = XLSX.readxlsx(filename)
        sheet = data["combined"]
        header = string.(vec(sheet[1, :]))
        rows = sheet[2:end, :]
        rows[.!ismissing.(rows)] = string.(rows[.!ismissing.(rows)])
        rows = convert.(Union{String, Missing}, rows)
        df_raw = DataFrame(rows, header)
    else
        df_raw = CSV.read(filename, DataFrame, stringtype = String)
        header = names(df_raw)
    end

    if any(ismissing.(df_raw[!, 1]))
        @warn "First column contains missing values"
    else
        df_raw[!, 1] = string.(df_raw[!, 1])
    end

    for col in header
        ET = eltype(df_raw[!, col])
        if ET <: Union{Missing, <:Number}
            continue
        end
        if nonmissingtype(ET) <: AbstractString
            if ET <: AbstractString
                if any(.!isnothing.(match.(r"^ *$", df_raw[!, col])))
                    @info "Converting $col to a missing union"
                    df_raw[!, col] = allowmissing(df_raw[!, col])
                    df_raw[.!isnothing.(match.(r"^ *$", df_raw[!, col])), col] .= missing
                end
            else
                for i in 1:nrow(df_raw)
                    val = df_raw[i, col]
                    if !ismissing(val) && !isnothing(match(r"^ *$", val))
                        df_raw[i, col] = missing
                    end
                end
            end
            if any(occursin.(r" [cm]?m$", skipmissing(df_raw[!, col])))
                df_raw[!, col] = uparse_missing.(df_raw[!, col])
            elseif !any(occursin.(r"[^0-9 ]", skipmissing(df_raw[!, col])))
                df_raw[!, col] = parse_missing.(Int, df_raw[!, col])
            elseif !any(occursin.(r"[^0-9. ]", skipmissing(df_raw[!, col])))
                df_raw[!, col] = parse_missing.(Float64, df_raw[!, col])
            elseif !any(occursin.(r"[^0-9, ]", skipmissing(df_raw[!, col])))
                df_raw[!, col] = parse_missing.(Vector{Int}, df_raw[!, col])
            elseif !any(occursin.(r"[^0-9., ]", skipmissing(df_raw[!, col])))
                df_raw[!, col] = parse_missing.(Vector{Float64}, df_raw[!, col])
            elseif !any(occursin.(r"[^0-9-, ]", skipmissing(df_raw[!, col])))
                df_raw[!, col] = parse_missing_range.(Int, df_raw[!, col])
            elseif !any(occursin.(r"[^0-9-., ]", skipmissing(df_raw[!, col])))
                df_raw[!, col] = parse_missing_range.(Float64, df_raw[!, col])
            elseif any(occursin.("2n", skipmissing(df_raw[!, col])))
                df_raw[!, col] = parse_missing_2n.(df_raw[!, col])
            elseif any(occursin.(",", skipmissing(df_raw[!, col])))
                df_raw[!, col] = parse_missing.(Vector{String}, df_raw[!, col])
            end
        end
    end

    return clean_traits(df_raw)
end

function read_plant_tree(path, synonyms, df)
    id = Dict{String, Int}()
    for (i, name) in enumerate(df.binomial)
        id[name] = i
    end

    tree = open(parsenewick, path)

    tips = getleafnames(tree)
    species_tips = replace.(tips, "_" => " ")
    for (old, new) in zip(tips, species_tips)
        renamenode!(tree, old, new)
    end

    renamed = 0
    tipset = Set(species_tips)
    duplicates = String[]
    for species in setdiff(species_tips, df.binomial)
        if species ∈ keys(synonyms)
            if synonyms[species] ∈ tipset
                push!(duplicates, species)
            else
                renamed += 1
                renamenode!(tree, species, synonyms[species])
            end
        end
    end
    droptips!(tree, duplicates)
    @info "Renamed $renamed species using synonyms, removed $(length(duplicates)) synonym duplicates"

    tips = getleafnames(tree)
    tipset = Set(tips)
    dropped_species = setdiff(tips, df.binomial)
    duplicates = String[]
    oddities = setdiff(getleafnames(tree), df.binomial)
    renamed = 0
    for species in oddities
        nn = split(species)
        if length(nn) == 3
            if (nn[1] * " " * nn[2] * "-" * nn[3]) ∈ df.binomial
                renamed += 1
                renamenode!(tree, species, nn[1] * " " * nn[2] * "-" * nn[3])
                continue
            elseif (nn[1] * " " * nn[2] * "-" * nn[3]) ∈ keys(synonyms)
                renamed += 1
                renamenode!(tree, species,
                            synonyms[nn[1] * " " * nn[2] * "-" * nn[3]])
                continue
            end
        end
        if length(nn) > 2
            if (nn[1] * " " * nn[2]) ∈ df.binomial
                if (nn[1] * " " * nn[2]) ∉ tipset
                    renamed += 1
                    renamenode!(tree, species, nn[1] * " " * nn[2])
                else
                    push!(duplicates, species)
                end
            elseif (nn[1] * " " * nn[2]) ∈ keys(synonyms)
                if synonyms[nn[1] * " " * nn[2]] ∉ tipset
                    renamed += 1
                    renamenode!(tree, species, synonyms[nn[1] * " " * nn[2]])
                else
                    push!(duplicates, species)
                end
            end
        end
    end
    droptips!(tree, duplicates)
    @info "Renamed $renamed abbreviated species, removed $(length(duplicates)) abbreviated synonym duplicates"

    tips = getleafnames(tree)
    dropped_species = setdiff(tips, df.binomial)
    @info dropped_species[1:min(10, length(dropped_species))]
    if !isempty(dropped_species)
        droptips!(tree, dropped_species)
        @warn "Dropped $(length(dropped_species)) species not in WFO"
    end

    ranks = ["genus", "family", "order", "phylum", "kingdom"]
    taxa = Set{String}((df[!, :] |> filter(row -> row.rank ∈ ranks)).binomial)
    @info "$(length(taxa)) taxa"

    setdiff!(taxa, Set(getnodenames(tree)))
    @info "$(length(taxa)) taxa not in parent tree"

    for node in traversal(tree, postorder)
        nd = getnodedata(tree, node)
        if isleaf(tree, node)
            ln = getnodename(tree, node)
            if !haskey(id, ln)
                @error "Could not find $ln"
            end
            row = df[id[ln], :]
            for rank in ranks
                nd[rank] = row[rank]
            end
        else
            children = getnodedata.(Ref(tree), getchildren(tree, node))
            for rank in ranks
                rvals = unique(get.(children, rank, 1:length(children)))
                if length(rvals) == 1
                    nd[rank] = rvals[1]
                end
            end
        end
    end

    for node in traversal(tree, preorder)
        if !isleaf(tree, node)
            nd = getnodedata(tree, node)
            for rank in ranks
                if haskey(nd, rank)
                    rval = nd[rank]
                    current = getnodename(tree, node)
                    if rval == current
                        break
                    else
                        if current ∉ taxa
                            if rval ∈ taxa
                                renamenode!(tree, node, rval)
                                delete!(taxa, rval)
                                break
                            end
                        else
                            @warn "Node $current has a real name, but should be $rval?"
                            break
                        end
                    end
                end
            end
        end
    end

    @info "Finally $(length(taxa)) taxa missing"

    n = length(taxa ∩ Set(df.genus))
    if n > 0
        @info "Missing $n of $(length(Set(df.genus))) genera"
    end
    n = length(taxa ∩ Set(df.family))
    if n > 0
        @info "Missing $n of $(length(Set(df.family))) families"
    end
    n = length(taxa ∩ Set(df.order))
    if n > 0
        @info "Missing $n of $(length(Set(df.order))) order"
    end
    n = length(taxa ∩ Set(df.phylum))
    if n > 0
        @info "Missing $n of $(length(Set(df.phylum))) phylum"
    end
    n = length(taxa ∩ Set(df.kingdom))
    if n > 0
        @info "Missing $n of $(length(Set(df.kingdom))) kingdom"
    end

    return tree
end

function expand_df_to_tree!(df, tree)
    set_spp = Set(df.species)
    for sp in getleafnames(tree)
        if sp ∉ set_spp
            push!(df, [sp; fill(missing, ncol(df) - 1)])
        end
    end
end

function expand_tree!(tree, df;
                      ranks = ["genus", "family", "order", "phylum", "kingdom"],
                      species = missing)
    leaves = Set(getleafnames(tree))
    num = length(leaves)
    height = getheight(tree, first(leaves))
    for row in eachrow(df)
        if row.rank == "species" && row.binomial ∉ leaves &&
           (ismissing(species) || row.binomial ∈ species)
            for rank in ranks
                if hasnode(tree, row[rank])
                    parent = getnode(tree, row[rank])
                    child = createnode!(tree, row.binomial)
                    for r in ["genus", "family", "order", "phylum", "kingdom"]
                        setnodedata!(tree, child, r, row[r])
                    end
                    createbranch!(tree, parent, child,
                                  height - getheight(tree, parent))
                    push!(leaves, row.binomial)
                    break
                end
            end
        end
    end
    @info "Added $(length(leaves) - num) new species, now $(length(leaves)) leaves"

    return tree
end

function synonymise(sheet, desynonym, sp_to_fam)
    cols = vec(string.(sheet[1, :]))
    species = [replace(sp, r"^ *([^ ].*[^ *])[ *]*$" => s"\1")
               for sp in vec(string.(sheet[2:end, 2]))]
    species = [get(desynonym, sp, sp) for sp in species]
    families = [get(sp_to_fam, sp, missing) for sp in species]
    if any(ismissing.(families))
        @warn "Missing $(length(species[ismissing.(families)])) records"
    end
    auth = vec(string.(sheet[2:end, 3]))
    ord = sortperm(species)

    return (vv = [families[ord], species[ord], auth[ord]], cols = cols)
end

function resynonym(species, synonyms)
    cols = ["Synonym", "Accepted name"]
    syn = String[]
    spp = String[]
    for sp in species
        for s in get(synonyms, sp, [])
            push!(syn, s)
            push!(spp, sp)
        end
    end
    return (vv = [syn, spp], cols = cols)
end

export read_trait_data, read_plant_tree, expand_tree!, expand_df_to_tree!,
       synonymise, resynonym

end
