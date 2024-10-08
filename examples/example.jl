using XLSX, DataFrames, CSV, JSON
using TraitExtraction
using Phylo

plantlist = JSON.parsefile("/Users/richardr/Downloads/World Flora Online/plant_list_2024-06.json");
# plantlist = JSON.parsefile("/Users/richardr/Downloads/World Flora Online/plant_list_2023-12.json");

df_dicot = read_trait_data("examples/dicots-6000.bhl.traits.xlsx");
df_monocot = read_trait_data("examples/angiosperm-monocot-1500.bhl.traits.xlsx");
df_0 = read_trait_data("../../plant-traits/output/raw/ocr-traits/angiosperm-wfo-0.bhl.traits.xlsx");

# Add columns to all dataframe so they have the same columns
cols = Set(names(df_dicot) ∪ names(df_monocot) ∪ names(df_0))
col_types = Set{Pair{String, Type}}()
for col in cols
    union_type = Union{}
    for dfi in [df_dicot, df_monocot, df_0]
        if col ∈ names(dfi)
            union_type = Union{union_type, eltype(dfi[!, col])}
        end
    end
    push!(col_types, col => union_type)
end

for (col, t) in col_types
    if col ∉ names(df_monocot)
        df_monocot[!, col] = t[missing for _ in 1:nrow(df_monocot)]
    end
    if col ∉ names(df_dicot)
        df_dicot[!, col] = t[missing for _ in 1:nrow(df_dicot)]
    end
    if col ∉ names(df_0)
        df_0[!, col] = t[missing for _ in 1:nrow(df_0)]
    end
end

df_full = vcat(df_dicot, df_monocot, df_0);
sort!(df_full, "taxon")
unique!(df_full)

species_genus = Dict{String, String}()
genus_family = Dict{String, String}()
synonym_to_species = Dict{String, String}()
species_to_synonyms = Dict{String, Set{String}}()
binomial_id = Dict{String, Int}()
name_id = Dict{String, Int}()
dropped = String[]

df_wfo = DataFrame(rank = String[], binomial = String[], name = String[],
                   species = Union{Nothing, String}[], genus = Union{Nothing, String}[],
                   family = Union{Nothing, String}[], order = Union{Nothing, String}[],
                   phylum = Union{Nothing, String}[], kingdom = Union{Nothing, String}[])
                   
for plant in plantlist
    if haskey(plant, "role_s")
        if plant["role_s"] == "accepted"
            row = (rank = plant["rank_s"], binomial = plant["full_name_string_alpha_s"],
                   name = plant["full_name_string_plain_s"],
                   species = get(plant, "placed_in_species_s", nothing),
                   genus = get(plant, "placed_in_genus_s", nothing),
                   family = get(plant, "placed_in_family_s", nothing),
                   order = get(plant, "placed_in_order_s", nothing),
                   phylum = get(plant, "placed_in_phylum_s", nothing),
                   kingdom = get(plant, "placed_in_kingdom_s", nothing))
            if row.rank != "code"
                push!(df_wfo, row)
                binomial_id[row.binomial] = nrow(df_wfo)
                name_id[row.name] = nrow(df_wfo)
                if row.rank == "species"
                    species_genus[row.binomial] = row.genus
                    species_to_synonyms[row.binomial] = Set{String}()
                elseif plant["rank_s"] == "genus"
                    genus_family[row.genus] = row.family
                end
            end
        end
    end
end

counts = Dict{String, Int}()
for rank in df_wfo.rank
    if rank ∈ ["species", "genus", "family", "order", "phylum", "kingdom"]
        counts[rank] = get(counts, rank, 0) + 1
    end
end
counts

non_species = 0
for plant in plantlist
    if haskey(plant, "role_s")
        if plant["role_s"] == "synonym"
            if plant["rank_s"] == "species"
                accepted = plant["accepted_full_name_string_plain_s"]
                name = plant["full_name_string_alpha_s"]
                if haskey(binomial_id, accepted)
                    row = df_wfo[binomial_id[accepted], :]
                    if row.rank == "species"
                        if name != accepted
                            synonym_to_species[name] = accepted
                            push!(species_to_synonyms[accepted], name)
                        end
                    elseif row.rank == "subspecies"
                        label = row.genus * " " * row.species
                        if name != label
                            synonym_to_species[name] = label
                            push!(species_to_synonyms[label], name)
                        end
                    else
                        non_species += 1
                    end
                elseif haskey(name_id, accepted)
                    row = df_wfo[name_id[accepted], :]
                    if row.rank == "species"
                        synonym_to_species[name] = row.binomial
                        push!(species_to_synonyms[row.binomial], name)
                    elseif row.rank == "subspecies" || row.rank == "variety" || row.rank == "form" || row.rank == "subvariety"
                        synonym_to_species[name] = row.genus * " " * row.species
                        push!(species_to_synonyms[row.genus * " " * row.species], name)
                    else
                        non_species += 1
                    end
                else
                    push!(dropped, accepted)
                end
            end
        end
    end
end

df = copy(df_full)
for col in names(df)
    if eltype(df[!, col]) ∉ [Union{String, Missing}, String]
        df[!, col] = [ismissing(x) ? missing : (x isa Int ? x : string(x)) for x in df[!, col]]
    end
end

for i in 1:nrow(df)
    if df[i, "taxon"] ∈ keys(synonym_to_species)
        df[i, "taxon"] = synonym_to_species[df[i, "taxon"]]
    end
end

sort!(df, "taxon")
spp = df.taxon
duplicates = spp[1:end-1] .== spp[2:end]
push!(duplicates, false)
delete!(df, duplicates)

# create a Dict of the unique non-missing values in each column with the column name as the key
trait_data = Dict{String, Union{Set{String}, Set{Int}}}()

for col in names(df)
    vals = unique(df[!, col])
    trait_data[col] = Set(skipmissing(vals))
end

max_len = maximum(length.(values(trait_data)))

# Drop any elements from trait_data whose entries have no "," in them
trait_single = Dict(k => v for (k, v) in trait_data if eltype(v) == Int || !any(occursin(",", x) for x in v))
trait_dict = Dict(k => v for (k, v) in trait_data if eltype(v) == String && any(occursin(",", x) for x in v))
trait_colour = Dict(k => v for (k, v) in trait_dict if eltype(v) == String && any(occursin("yellow", x) for x in v))
trait_dict = Dict(k => v for (k, v) in trait_dict if eltype(v) == String && !any(occursin("yellow", x) for x in v))


ks_multi = [keys(trait_dict)...]
col_multi = Vector{Vector{Union{Missing,String}}}()

for k in ks_multi
    v = trait_dict[k]
    vals = Union{Missing, String}[missing for _ in 1:max_len]
    vals[1:length(v)] .= sort([v...])
    push!(col_multi, vals)
end

ks_multi = ks_multi[sortperm(length.(values(trait_dict)))]
col_multi = col_multi[sortperm(length.(values(trait_dict)))]

ks_colour = [keys(trait_colour)...]
col_colour = Vector{Vector{Union{Missing,String}}}()

for k in ks_colour
    v = trait_colour[k]
    vals = Union{Missing, String}[missing for _ in 1:max_len]
    vals[1:length(v)] .= sort([v...])
    push!(col_colour, vals)
end

ks_colour = ks_colour[sortperm(length.(values(trait_colour)))]
col_colour = col_colour[sortperm(length.(values(trait_colour)))]

ks_single = [keys(trait_single)...]
col_single = Vector{Vector{Union{Missing,String, Int}}}()

for k in ks_single
    v = trait_single[k]
    vals = Union{Missing, eltype(v)}[missing for _ in 1:max_len]
    vals[1:length(v)] .= sort([v...])
    push!(col_single, vals)
end

ks_single = ks_single[sortperm(length.(values(trait_single)))]
col_single = col_single[sortperm(length.(values(trait_single)))]

# Write the data to an excel file
XLSX.writetable("examples/traits.xlsx",
                [("Multi", col_multi, ks_multi),
                 ("Colour", col_colour, ks_colour),
                 ("Single", col_single, ks_single)],
                 overwrite = true)

full_tree = read_plant_tree("examples/Qian2016.tree", synonym_to_species, df_wfo);
expand_tree!(full_tree, df_wfo);

tree = read_plant_tree("examples/Qian2016.tree", synonym_to_species, df_wfo);
expand_tree!(tree, df_wfo, ranks = ["genus", "family"]);

real_tree = read_plant_tree("examples/Qian2016.tree", synonym_to_species, df_wfo);

df_tree = copy(df)
for sp in df_tree.taxon
    if !haskey(species_genus, sp)
        println(sp ∈ [x["full_name_string_alpha_s"] for x in plantlist] ?
                "$sp, in WFO ($([x["role_s"] for x in plantlist if sp == x["full_name_string_alpha_s"]][1]))" : "$sp, not in WFO")
        if sp ∈ getleafnames(tree)
            @warn "Species $sp in tree"
        else
            @info "Species $sp not in tree"
            delete!(df_tree, findall(df_tree.taxon .== sp))
        end
    end
end

dropped_species = setdiff(getleafnames(real_tree), df_tree.taxon)
droptips!(real_tree, dropped_species)
delete!(df_tree, df_tree.taxon .∉ Ref(getleafnames(real_tree)))
sort!(df_tree, "taxon")
spp = df_tree.taxon
duplicates = spp[1:end-1] .== spp[2:end]
push!(duplicates, false)
delete!(df_tree, duplicates)

using RCall
@rput real_tree
val = df_tree[:, "life form"]
@rput val
R"""
library(ape)
output <- ace(val, real_tree, "discrete", model = "SYM")
"""


species_wfo = df_wfo |> filter(row -> row.rank == "species")
select!(species_wfo, Not("rank"))

# Find duplicated binomials
spp = species_wfo.binomial
duplicates = spp[1:end-1] .== spp[2:end]
push!(duplicates, false)
delete!(species_wfo, duplicates)
select!(species_wfo, Not("name", "species"))
CSV.write("examples/wfo_species.csv", species_wfo)
# sp0 = CSV.read("../../plant-traits/input/wfo-species/wfo-0.csv", DataFrame)

data = XLSX.readxlsx("examples/Full PUP II sample.xlsx")
bryo_sheet = data["Bryophyte species"]
pter_sheet = data["Pteridophyte species"]
monocot_sheet = data["Monocot species"]
dicot_sheet = data["Dicot species"]

bryo_sp = vec(string.(bryo_sheet[2:end, 2]))
pter_sp = vec(string.(pter_sheet[2:end, 2]))
monocot_sp = vec(string.(monocot_sheet[2:end, 2]))
dicot_sp = vec(string.(dicot_sheet[2:end, 2]))

bryo_new = Dict(sp => synonym_to_species[sp] for sp in bryo_sp if sp in keys(synonym_to_species))
pter_new = Dict(sp => synonym_to_species[sp] for sp in pter_sp if sp in keys(synonym_to_species))
monocot_new = Dict(sp => synonym_to_species[sp] for sp in monocot_sp if sp in keys(synonym_to_species))
dicot_new = Dict(sp => synonym_to_species[sp] for sp in dicot_sp if sp in keys(synonym_to_species))
