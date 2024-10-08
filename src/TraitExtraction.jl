module TraitExtraction

using XLSX, DataFrames, Unitful, Statistics, Phylo

function uparse_missing(::Missing)
    return missing
end

function uparse_missing(val::String)
    return uparse(replace(val, r" " => ""))
end

function parse_missing(::Missing)
    return missing
end

function parse_missing(val::String)
    if occursin("-", val)
        nums = parse.(Float64, split(val, "-"))
        if length(nums) == 2
            return round(Int, mean(nums))
        else
            @warn "Could not parse $val"
            return missing
        end
    else
        return round(Int, parse(Float64, val))
    end
end

function clean_traits(df_raw)
    df = copy(df_raw)

    col = df[!, "flower architecture"]
    col[ismissing.(col) .| (col .== "catkin")] .= missing
    col[ismissing.(col) .| (col .== "pentamerous")] .= missing

    col = df[!, "fruit structure"]
    col[ismissing.(col) .| (col .== "dry, fleshy")] .= missing

    col = df[!, "petal fusion"]
    col[.!ismissing.(col) .& (col .== "free, fused")] .= "fused"

    col = df[!, "spinescence"]
    col[ismissing.(col) .| (col .== "armed, spinescent, spinulate, unarmed")] .= missing

    col = df[!, "fruit dehiscence"]
    col[.!ismissing.(col) .& (col .== "dehiscent, indehiscent, valvate")] .= "valvate"
    
    col = df[!, "habitat"]
    col[.!ismissing.(col) .& (col .== "aquatic, aquatic/hygrophilous")] .= "hygrophilous"
    col[.!ismissing.(col) .& (col .== "epiphyte, terrestrial")] .= "epiphyte"
  
    col = df[!, "succulence"]
    col[.!ismissing.(col) .& (col .== "flexhy")] .= "fleshy"
 
    col = df[!, "flower sex"]
    col[.!ismissing.(col) .& (col .== "bisexual, pistillate, staminate")] .= "unisexual"
    col[.!ismissing.(col) .& (col .== "bisexual, pistillate, unisexual")] .= "bissexual, pistillate"
    col[.!ismissing.(col) .& (col .== "pistillate, staminate, unisexual")] .= "unisexual"
    col[.!ismissing.(col) .& (col .== "staminate, unisexual")] .= "staminate"
    col[.!ismissing.(col) .& (col .== "pistillate, unisexual")] .= "unisexual"
 
    col = df[!, "leaf arrangement"]
    col[.!ismissing.(col) .& (col .== "bipinnate, lobed, pinnate")] .= "bipinnate"
    col[.!ismissing.(col) .& (col .== "bipinnate, pinnate")] .= "bipinnate"
    col[.!ismissing.(col) .& (col .== "digitate, pinnate")] .= "digitate"
    col[.!ismissing.(col) .& (col .== "lobed, pinnate")] .= "pinnate"
    col[.!ismissing.(col) .& (col .== "lobed")] .= "simple"
    col[.!ismissing.(col) .& (col .== "lobed, simple")] .= "simple"
    col[ismissing.(col) .| (col .== "lobed, pinnate, simple")] .= missing
    col[ismissing.(col) .| (col .== "simple, trifoliolate")] .= missing

    col = df[!, "leaf position"]
    col[ismissing.(col) .| (col .== "alternate, opposite")] .= missing
    col[ismissing.(col) .| (col .== "alternate, fascicled, opposite")] .= missing
    col[.!ismissing.(col) .& (col .== "decussate, opposite")] .= "decussate"
    col[.!ismissing.(col) .& (col .== "alternate, spirally")] .= "spirally"
    col[.!ismissing.(col) .& (col .== "fascicled")] .= "fasciculate"
    col[.!ismissing.(col) .& (col .== "verticillate, whorled")] .= "whorled"
    col[ismissing.(col) .| (col .== "alternate, fasciculate")] .= missing

    col = df[!, "fruit shape"]
    col[.!ismissing.(col) .& (col .== "cylindrical, fusiform")] .= "fusiform"
    col[.!ismissing.(col) .& (col .== "cylindrical, ovoid")] .= "ovoid"
    col[.!ismissing.(col) .& (col .== "turbinate")] .= "ovoid"

    col = df[!, "inflorescence arrangement"]
    col[.!ismissing.(col) .& (col .== "clustered, corymb, raceme")] .= "clustered, corymb"
    col[.!ismissing.(col) .& (col .== "clustered, panicle, raceme")] .= "clustered, panicle"
    col[.!ismissing.(col) .& (col .== "raceme, thyrse")] .= "thyrse"
    col[.!ismissing.(col) .& (col .== "raceme, umbel")] .= "umbel"
    col[.!ismissing.(col) .& (col .== "panicle, thyrse")] .= "thyrse"
    col[.!ismissing.(col) .& (col .== "corymb, panicle, raceme")] .= "corymb, panicle"

    col = df[!, "habit"]
    for i in 1:nrow(df)
        if ismissing(col[i])
            continue
        end
        if occursin(r"^climber", col[i])
            col[i] = "climber"
        elseif occursin(r"^character", col[i])
            col[i] = "climber"
        end
    end
    col[.!ismissing.(col) .& (col .== "erect leafy, herb")] .= "erect leafy"
    col[.!ismissing.(col) .& (col .== "erect leafy, herb, shrub")] .= "erect leafy"
    col[.!ismissing.(col) .& (col .== "herb, prostrate")] .= "prostrate"
    col[.!ismissing.(col) .& (col .== "herb, prostrate, shrub")] .= "prostrate"
    col[.!ismissing.(col) .& (col .== "herb, tree, tree/shrub")] .= "herb"
    col[.!ismissing.(col) .& (col .== "tree, tussock")] .= "tree"
    col[.!ismissing.(col) .& (col .== "shrub, tussock")] .= "shrub"
    col[.!ismissing.(col) .& (col .== "prostrate, shrub, tree/shrub")] .= "prostrate"
    col[.!ismissing.(col) .& (col .== "prostrate, tussock")] .= "prostrate"
    col[.!ismissing.(col) .& (col .== "herb, tussock")] .= "erect leafy"

    return df
end

function read_trait_data(filename::String)
    data = XLSX.readxlsx(filename)
    sheet = data["combined"]
    header = string.(vec(sheet[1, :]))
    rows = sheet[2:end, :]
    rows[.!ismissing.(rows)] = string.(rows[.!ismissing.(rows)])
    rows = convert.(Union{String, Missing}, rows)
    df_raw = DataFrame(rows, header)

    if any(ismissing.(df_raw[!, 1]))
        @warn "First column contains missing values"
    else
        df_raw[!, 1] = string.(df_raw[!, 1])
    end

    df = clean_traits(df_raw)

    for col in header
        if all(ismissing.(df[!, col]))
            select!(df, Not(col))
            @info "Dropping column $col"
            continue
        end
        if any(occursin.(r" [cm]?m$", skipmissing(df[!, col])))
            df[!, col] = uparse_missing.(df[!, col])
        elseif !any(occursin.(r"[^0-9-. ]", skipmissing(df[!, col])))
            df[!, col] = parse_missing.(df[!, col])
        else
            for i in 1:nrow(df)
                val = df[i, col]
                if ismissing(val)
                    continue
                end
                if occursin(r"2n", col)
                    val = replace(val, r"[,.] ?$" => "", r" ?2n ?= ?" => "", r" ?\| ?" => ", ", r" ?, ?" => ", ")
                    if occursin(r"[^0-9, ]", val) || !occursin(r"[^ ]", val)
                        @warn "Could not parse $(df[i, col])"
                        df[i, col] = missing
                    else
                        df[i, col] = join(string.(unique(sort(parse.(Int, split(val, ", "))))), ", ")
                    end
                elseif occursin(r" ?, ?", val)
                    df[i, col] = join(sort(split(val, r" ?, ?")), ", ")
                end
            end
        end
    end
    
    return df
end

function read_plant_tree(path, synonyms, df)
    id = Dict{String, Int}()
    for (i, name) in enumerate(df.binomial)
        id[name] = i
    end

    tree = open(parsenewick, path);

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
                renamenode!(tree, species, synonyms[nn[1] * " " * nn[2] * "-" * nn[3]])
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
    @info "$(length(taxa)) taxa not used"

    @info "$(nnodes(tree)) nodes"
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

    @info "Then $(length(taxa))"

    @info "$(length(taxa ∩ Set(df.genus))) of $(length(Set(df.genus))) genera"
    @info "$(length(taxa ∩ Set(df.family))) of $(length(Set(df.family))) families"
    @info "$(length(taxa ∩ Set(df.order))) of $(length(Set(df.order))) order"
    @info "$(length(taxa ∩ Set(df.phylum))) of $(length(Set(df.phylum))) phylum"
    @info "$(length(taxa ∩ Set(df.kingdom)))  of $(length(Set(df.kingdom))) kingdom"

    return tree
end

function expand_tree!(tree, df; ranks = ["genus", "family", "order", "phylum", "kingdom"])
    leaves = Set(getleafnames(tree))
    height = getheight(tree, first(leaves))
    for row in eachrow(df)
        if row.rank == "species" && row.binomial ∉ leaves
            for rank in ranks
                if hasnode(tree, row[rank])
                    parent = getnode(tree, row[rank])
                    child = createnode!(tree, row.binomial)
                    for r in ["genus", "family", "order", "phylum", "kingdom"]
                        setnodedata!(tree, child, r, row[r])
                    end
                    createbranch!(tree, parent, child, height - getheight(tree, parent))
                    push!(leaves, row.binomial)
                    break
                end
            end
        end
    end

    return tree
end

export read_trait_data, read_plant_tree, expand_tree!

end
