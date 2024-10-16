using XLSX, DataFrames, CSV, JSON
using TraitExtraction
using Phylo
using DataStructures
using PhyloNetworks, StatsModels

plantlist = JSON.parsefile("/Users/richardr/Downloads/World Flora Online/plant_list_2024-06.json");
# plantlist = JSON.parsefile("/Users/richardr/Downloads/World Flora Online/plant_list_2023-12.json");

df_all = read_trait_data("data/combined.csv");
rename!(df_all, :taxon => :species)

df_full = copy(df_all)
sort!(df_full, "species")
unique!(df_full)
counter(eltype.(eachcol(df_full)))

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
                   phylum = Union{Nothing, String}[], kingdom = Union{Nothing, String}[],
                   wfo = Union{Nothing, String}[])
                   
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
                   kingdom = get(plant, "placed_in_kingdom_s", nothing),
                   wfo = haskey(plant, "wfo_id_s") ?
                        "https://www.worldfloraonline.org/taxon/" * plant["wfo_id_s"] : nothing)
            if row.rank != "code" && !haskey(binomial_id, row.binomial)
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

counts = Dict(counter(df_wfo.rank))
ranks = collect(keys(counts))
for rank in ranks
    if rank ∉ ["species", "genus", "family", "order", "phylum", "kingdom"]
        delete!(counts, rank)
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
                        if name != accepted && !haskey(binomial_id, name)
                            synonym_to_species[name] = accepted
                            push!(species_to_synonyms[accepted], name)
                        end
                    elseif row.rank == "subspecies"
                        label = row.genus * " " * row.species
                        if name != label && !haskey(binomial_id, name)
                            synonym_to_species[name] = label
                            push!(species_to_synonyms[label], name)
                        end
                    else
                        non_species += 1
                    end
                elseif haskey(name_id, accepted)
                    row = df_wfo[name_id[accepted], :]
                    if row.rank == "species"
                        if name != row.binomial && !haskey(binomial_id, name)
                            synonym_to_species[name] = row.binomial
                            push!(species_to_synonyms[row.binomial], name)
                        end
                    elseif row.rank == "subspecies" || row.rank == "variety" || row.rank == "form" || row.rank == "subvariety"
                        label = row.genus * " " * row.species
                        if name != label && !haskey(binomial_id, name)
                            synonym_to_species[name] = label
                            push!(species_to_synonyms[label], name)
                        end
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

df_full.species = get.(Ref(synonym_to_species), df_full.species, df_full.species)
sort!(df_full, "species")
spp = df_full.species;
duplicates = spp[1:end-1] .== spp[2:end];
push!(duplicates, false);
delete!(df_full, duplicates)
@info "Dropped $(sum(duplicates)) duplicates"

# Create tree here
df = copy(df_full)
for col in names(df)
    if !(nonmissingtype(eltype(df[!, col])) <: AbstractString)
        if nonmissingtype(eltype(df[!, col])) <: Vector{<:AbstractString}
            df[!, col] = [ismissing(val) ? missing : join(val, ",") for val in df[!, col]]
        else
            df[!, col] = [ismissing(val) ? missing : string(val) for val in df[!, col]]
        end
    end
end

# create a Dict of the unique non-missing values in each column with the column name as the key
trait_data = Dict{String, Union{Set{String}, Set{Int}}}()

for col in names(df)
    if col ∉ ["wfo", "species"]
        if !(nonmissingtype(eltype(df_full[!, col])) <: Number ||
            (nonmissingtype(eltype(df_full[!, col])) <: AbstractVector &&
             eltype(nonmissingtype(eltype(df_full[!, col]))) <: Number))
            @info col
            cnt = Dict(counter(df[!, col]))
            for (k, v) in cnt
                if v < 50
                    if length(unique(df_full[!, col])) < 20
                        @info "Deleting $col: $k with $v occurrences"
                    end
                    delete!(cnt, k)
                    TraitExtraction.removelabel!(df[!, col], k)
                end
            end
            if length(keys(cnt)) ≤ 2
                @info "Deleting $col: all missing or only 1 value"
                select!(df, Not(col))
            else
                if length(unique(df[!, col])) < 20
                    println(counter(df[!, col]))
                else
                    println("$(length(unique(df[!, col]))) unique values")
                end
                trait_data[col] = Set(skipmissing(unique(df[!, col])))
            end
        end
    end
end

df.wfo = allowmissing(df.wfo);
df.detected = allowmissing(df.detected);
sort!(df, "species")

real_tree = read_plant_tree("data/Qian2016.tree", synonym_to_species, df_wfo);
expand_df_to_tree!(df, real_tree);
open("examples/real_tree.tree", "w") do io
    write(io, real_tree)
end

max_len = maximum(length.(values(trait_data)))

#=
full_tree = read_plant_tree("data/Qian2016.tree", synonym_to_species, df_wfo);
expand_tree!(full_tree, df_wfo);

tree = read_plant_tree("data/Qian2016.tree", synonym_to_species, df_wfo);
expand_tree!(tree, df_wfo, ranks = ["genus", "family"]);
=#

df_tree = copy(df)
for sp in df_tree.species
    if !haskey(species_genus, sp)
        println(sp ∈ [x["full_name_string_alpha_s"] for x in plantlist] ?
                "$sp, in WFO ($([x["role_s"] for x in plantlist if sp == x["full_name_string_alpha_s"]][1]))" : "$sp, not in WFO")
        if sp ∈ getleafnames(real_tree)
            @warn "Species $sp in tree"
        else
            @info "Species $sp not in tree"
            delete!(df_tree, findall(df_tree.species .== sp))
        end
    end
end

expand_tree!(real_tree, df_wfo, ranks = ["genus", "family"], species = Set(df_tree.species))
# expand_tree!(real_tree, df_wfo, ranks = ["genus", "family"], species = df_tree.species)
expand_df_to_tree!(df_tree, real_tree);
sort!(df_tree, "species")
spp = df_tree.species
duplicates = spp[1:end-1] .== spp[2:end]
push!(duplicates, false)
delete!(df_tree, duplicates)

dropped_species = setdiff(getleafnames(real_tree), df_tree.species)
if !isempty(dropped_species)
    println("Dropped $(length(dropped_species)) tips from tree")
    droptips!(real_tree, dropped_species)
end
delete!(df_tree, df_tree.species .∉ Ref(Set(getleafnames(real_tree))))
nleaves(real_tree)

col = "leaf min. length [cm]";
df_pn = df_tree[!, ["species", col]];
df_pn[!, "tipNames"] = replace.(df_tree.species, " " => "_");
# rename!(df_pn, col => "trait")
df_pn[!, "trait"] = TraitExtraction.uparse_missing.(df_pn[!, col]);
df_pn[!, "trait"] ./= oneunit(nonmissingtype(eltype(df_full[!, col])));
sum(skipmissing(df_pn[!, "trait"]) .< 200)
for i in eachindex(df_pn[!, "trait"])
    if !ismissing(df_pn.trait[i]) && (df_pn.trait[i] < 0.1 || df_pn.trait[i] > 200)
        df_pn.trait[i] = missing
    end
    if !ismissing(df_pn.trait[i])
        df_pn.trait[i] = log(df_pn.trait[i])
    end
end

cnt = counter(getindex.(split.(df_full.species), 1));
genera = collect(keys(cnt))[sortperm(collect(values(cnt)))][6650:end];
species = df_wfo.binomial[(df_wfo.genus .∈ Ref(Set(genera))) .& (df_wfo.rank .== "species")]
expand_tree!(real_tree, df_wfo, ranks = ["genus"], species = Set(species))

delete!(df_pn, findall(ismissing.(df_pn[!, "trait"]) .& (df_pn.species .∉ Ref(Set(species)))))
keeptips!(real_tree, df_pn.species)
select!(df_pn, Not(["species", col]))

open("examples/tree.tree", "w") do io
    write(io, real_tree)
end

# sed 's/"//g' tree.tree > minus.tree
# sed 's/ /_/g' minus.tree > final.tree

pntree = readTopology("examples/final.tree");

fitBM = phylolm(@formula(trait ~ 1), df_pn, pntree);
ancStates = ancestralStateReconstruction(fitBM); # Should produce a warning, as variance is unknown.

mnames = Set(df_pn.tipNames[ismissing.(df_pn.trait)]);
ex = expectations(ancStates);
ex.condExpectation = exp.(ex.condExpectation) .* oneunit(nonmissingtype(eltype(df_full[!, col])))
ex.lower = exp.(predint(ancStates)[:,1]) .* oneunit(nonmissingtype(eltype(df_full[!, col])))
ex.upper = exp.(predint(ancStates)[:,2]) .* oneunit(nonmissingtype(eltype(df_full[!, col])))

ex[ex.nodeNumber .∈ Ref(Set(findall(pntree.names .∈ Ref(mnames)))), :]

pip = predintPlot(ancStates);
pip[pip.nodeNumber .∈ Ref(Set(findall(pntree.names .∈ Ref(mnames)))), :]

using RCall
R"""
library(ape)
library(expm)
"""
@rput real_tree
R"""
tips <- real_tree$tip.label
"""
@rget tips
ord = sortperm(tips)

states_df = df_tree[sortperm(ord), ["species", "life form"]]
all(states_df.species .== tips)
@rput states_df

R"""
output <- ace(states_df[["life form"]], real_tree, "discrete", model = "SYM")
rates <- output$rates[output$index.matrix]
rates <- matrix(rates, dim(output$index.matrix)[1], dim(output$index.matrix)[2])
for (i in 1:dim(rates)[1]) {
    rates[i, i] = 0
    rates[i, i] = -sum(rates[i,])
}
for (n in 1:dim(real_tree$edge)[1]) {
    results <- states_df[["life form"]]
    if (real_tree$edge[n, 2] <= length(states) && is.na(states[n])) {
        tprobs <- expm(rates * real_tree$edge.length[n])
        node <- real_tree$edge[n, 1] - min(real_tree$edge[,1]) + 1
        probs <- tprobs %*% t(t(output$lik.anc[node,]))
        idx <- which.max(probs)
        if (probs[idx] > 0.5) {
            results[real_tree$edge[n, 2]] <- colnames(output$lik.anc)[idx]
        }
    }
}
"""
@rget results
counter(states_df[!, "life form"])
counter(results)
@rget rates
df_data = copy(df_tree)
rows = findall(ismissing.(df_data[!, 2]))
@info 2, length(rows)
for col in 3:ncol(df_data)
    rows = rows[findall(ismissing.(df_data[rows, col]))]
    @info col, length(rows)
end
df_data = df_data[setdiff(1:nrow(df_data), rows), :]
reduced_tree = read_plant_tree("data/Qian2016.tree", synonym_to_species, df_wfo);
dropped_species = setdiff(getleafnames(reduced_tree), df_data.species)
if !isempty(dropped_species)
    println("Dropped $length(dropped_species) tips from tree")
    droptips!(reduced_tree, dropped_species)
end
vals = TraitExtraction.uparse_missing.(df_data[:, "plant min. height [m]"]) ./ oneunit(eltype(df_full[!, "plant min. height [m]"]))
@rput vals
@rput reduced_tree

R"""
library(ape)
library(expm)
output_c <- ace(vals, reduced_tree, "continuous")
cb <- corBrownian(1, reduced_tree, form = ~1)
output_c <- ace(vals, reduced_tree, method = "GLS", corStruct = cb)
rates <- output$rates[output$index.matrix]
rates <- matrix(rates, dim(output$index.matrix)[1], dim(output$index.matrix)[2])
for (i in 1:dim(rates)[1]) {
    rates[i, i] = 0
    rates[i, i] = -sum(rates[i,])
}
for (n in 1:dim(reduced_tree$edge)[1]) {
    results <- states
    if (reduced_tree$edge[n, 2] <= length(states) && is.na(states[n])) {
        tprobs <- expm(rates * reduced_tree$edge.length[n])
        node <- reduced_tree$edge[n, 1] - min(reduced_tree$edge[,1]) + 1
        probs <- tprobs %*% t(t(output$lik.anc[node,]))
        idx <- which.max(probs)
        if (probs[idx] > 0.5) {
            results[reduced_tree$edge[n, 2]] <- colnames(output$lik.anc)[idx]
        } else {
            print(c(reduced_tree$edge[n, 2], idx, probs[idx]))
        }
    }
}
"""
@rget results

species_wfo = df_wfo |> filter(row -> row.rank == "species")
select!(species_wfo, Not("rank"))

# Find duplicated binomials
spp = species_wfo.binomial
duplicates = spp[1:end-1] .== spp[2:end]
push!(duplicates, false)
delete!(species_wfo, duplicates)
select!(species_wfo, Not("name", "species"))
CSV.write("examples/wfo_species.csv", species_wfo)

species_to_family = Dict{String, String}()
for sp in keys(species_genus)
    species_to_family[sp] = genus_family[species_genus[sp]]
end

data = XLSX.readxlsx("data/Full PUP II sample.xlsx");
bryo = synonymise(data["Bryophyte species"], synonym_to_species, species_to_family);
bryo_syn = resynonym(bryo.vv[2], species_to_synonyms);
pter = synonymise(data["Pteridophyte species"], synonym_to_species, species_to_family);
pter_syn = resynonym(pter.vv[2], species_to_synonyms);
monocot = synonymise(data["Monocot species"], synonym_to_species, species_to_family);
monocot_syn = resynonym(monocot.vv[2], species_to_synonyms);
dicot = synonymise(data["Dicot species"], synonym_to_species, species_to_family);
dicot_syn = resynonym(dicot.vv[2], species_to_synonyms);

XLSX.writetable("data/New PUP II.xlsx",
                [("Bryophyte species", bryo.vv, bryo.cols),
                 ("Bryophyte synonyms", bryo_syn.vv, bryo_syn.cols),
                 ("Pteridophyte species", pter.vv, pter.cols),
                 ("Pteridophyte synonyms", pter_syn.vv, pter_syn.cols),
                 ("Monocot species", monocot.vv, monocot.cols),
                 ("Monocot synonyms", monocot_syn.vv,monocot_syn.cols),
                 ("Dicot species", dicot.vv, dicot.cols),
                 ("Dicot synonyms", dicot_syn.vv, dicot_syn.cols)],
                 overwrite = true)

# Drop any elements from trait_data whose entries have no "," in them
trait_single = Dict(k => v for (k, v) in trait_data if eltype(v) <: Number || !any(occursin(",", x) for x in v))
trait_dict = Dict(k => v for (k, v) in trait_data if eltype(v) <: AbstractString && any(occursin(",", x) for x in v))
trait_colour = Dict(k => v for (k, v) in trait_dict if eltype(v) <: AbstractString && any(occursin("yellow", x) for x in v))
trait_dict = Dict(k => v for (k, v) in trait_dict if eltype(v) <: AbstractString && !any(occursin("yellow", x) for x in v))


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

climate = CSV.read("data/climate.csv", DataFrame)
for col in names(climate)[2:end]
    climate[!, col] = passmissing(x -> round(x, sigdigits = 4)).(climate[!, col])
end

# Add WFO refs into climate data and round to 4dp
climate.wfo = [haskey(binomial_id, row.species) ? df_wfo.wfo[binomial_id[row.species]] : missing for row in eachrow(climate)]
delete!(climate, findall(ismissing.(climate.wfo)))
CSV.write("examples/climate.csv", climate)
length(climate.species), length(df_full.species), length(climate.species ∩ df_full.species)
