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

function relabel!(col, old, new)
    col[.!ismissing.(col) .& (col .== old)] .= new
    return col
end

function clean_traits(df_raw)
    df = copy(df_raw)

    col = df[!, "flower architecture"]
    relabel!(col, "catkin", missing)
    relabel!(col, "pentamerous", missing)

    col = df[!, "fruit structure"]
    relabel!(col, "dry, fleshy", missing)

    col = df[!, "petal fusion"]
    relabel!(col, "free, fused", missing)
    relabel!(col, "free, fusion", missing)
    relabel!(col, "fused, fusion", "fused")

    col = df[!, "spinescence"]
    relabel!(col, "armed, spinescent, spinulate, unarmed", "spinescent, spinulate")
    relabel!(col, "armed, spinescent, spinulate", "spinescent, spinulate")
    relabel!(col, "armed, spinescent", "spinescent")
    relabel!(col, "armed, spinulate", "spinulate")
    relabel!(col, "armed, unarmed", missing)
    relabel!(col, "spinescent, unarmed", "spinescent")
    relabel!(col, "spinulate, unarmed", "spinulate")

    col = df[!, "fruit dehiscence"]
    relabel!(col, "dehiscent, indehiscent, valvate", "valvate")
    relabel!(col, "dehiscent, suture", "suture")
    
    col = df[!, "habitat"]
    relabel!(col, "aquatic, aquatic/hygrophilous", "hygrophilous")
    relabel!(col, "epiphyte, terrestrial", "epiphyte")
  
    col = df[!, "succulence"]
    relabel!(col, "flexhy", "fleshy")
 
    col = df[!, "flower sex"]
    relabel!(col, "bisexual, pistillate, staminate", "unisexual")
    relabel!(col, "bisexual, pistillate, unisexual", "bissexual, pistillate")
    relabel!(col, "pistillate, staminate, unisexual", "unisexual")
    relabel!(col, "staminate, unisexual", "staminate")
    relabel!(col, "pistillate, unisexual", "unisexual")
 
    col = df[!, "leaf arrangement"]
    relabel!(col, "bipinnate, lobed, pinnate", "bipinnate")
    relabel!(col, "bipinnate, pinnate", "bipinnate")
    relabel!(col, "digitate, pinnate", "digitate")
    relabel!(col, "lobed, pinnate", "pinnate")
    relabel!(col, "lobed", "simple")
    relabel!(col, "lobed, simple", "simple")
    relabel!(col, "lobed, pinnate, simple", missing)
    relabel!(col, "simple, trifoliolate", missing)

    col = df[!, "leaf position"]
    relabel!(col, "alternate, opposite", missing)
    relabel!(col, "alternate, fascicled, opposite", missing)
    relabel!(col, "decussate, opposite", "decussate")
    relabel!(col, "alternate, spirally", "spirally")
    relabel!(col, "fascicled", "fasciculate")
    relabel!(col, "verticillate, whorled", "whorled")
    relabel!(col, "alternate, fasciculate", missing)

    col = df[!, "fruit shape"]
    relabel!(col, "cylindrical, fusiform", "fusiform")
    relabel!(col, "cylindrical, ovoid", "ovoid")
    relabel!(col, "turbinate", "ovoid")

    col = df[!, "inflorescence arrangement"]
    relabel!(col, "clustered, corymb, raceme", "clustered, corymb")
    relabel!(col, "clustered, panicle, raceme", "clustered, panicle")
    relabel!(col, "clustered, panicle, raceme, thyrse", "clustered, thyrse")
    relabel!(col, "raceme, thyrse", "thyrse")
    relabel!(col, "raceme, umbel", "umbel")
    relabel!(col, "panicle, thyrse", "thyrse")
    relabel!(col, "corymb, panicle, raceme", "corymb, panicle")
    relabel!(col, "panicle, raceme, umbel", "panicle, umbel")

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
    relabel!(col, "erect leafy, herb", "erect leafy")
    relabel!(col, "erect leafy, herb, shrub", "erect leafy")
    relabel!(col, "herb, prostrate", "prostrate")
    relabel!(col, "herb, prostrate, shrub", "prostrate")
    relabel!(col, "herb, tree, tree/shrub", "herb")
    relabel!(col, "tree, tussock", "tree")
    relabel!(col, "shrub, tussock", "shrub")
    relabel!(col, "prostrate, shrub, tree/shrub", "prostrate")
    relabel!(col, "prostrate, tussock", "prostrate")
    relabel!(col, "herb, tussock", "erect leafy")

    for col in names(df)
        if all(ismissing.(df[!, col]))
            select!(df, Not(col))
            @info "Dropping column $col"
        end
    end

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

    for col in header
        if all(ismissing.(df_raw[!, col]))
            continue
        end
        if any(occursin.(r" [cm]?m$", skipmissing(df_raw[!, col])))
            df_raw[!, col] = uparse_missing.(df_raw[!, col])
        elseif !any(occursin.(r"[^0-9-. ]", skipmissing(df_raw[!, col])))
            df_raw[!, col] = parse_missing.(df_raw[!, col])
        else
            for i in 1:nrow(df_raw)
                val = df_raw[i, col]
                if ismissing(val)
                    continue
                end
                if occursin(r"2n", col)
                    val = replace(val, r"[,.] ?$" => "", r" ?2n ?= ?" => "", r" ?\| ?" => ", ", r" ?, ?" => ", ")
                    if occursin(r"[^0-9, ]", val) || !occursin(r"[^ ]", val)
                        @warn "Could not parse $(df_raw[i, col])"
                        df_raw[i, col] = missing
                    else
                        df_raw[i, col] = join(string.(unique(sort(parse.(Int, split(val, ", "))))), ", ")
                    end
                elseif occursin(r" ?, ?", val)
                    df_raw[i, col] = join(sort(split(val, r" ?, ?")), ", ")
                end
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

function expand_tree!(tree, df; ranks = ["genus", "family", "order", "phylum", "kingdom"])
    leaves = Set(getleafnames(tree))
    num = length(leaves)
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
    @info "Added $(length(leaves) - num) new species, now $(length(leaves)) leaves"

    return tree
end

export read_trait_data, read_plant_tree, expand_tree!

end
