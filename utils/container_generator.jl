#!/usr/bin/env julia
using SimpleContainerGenerator

mkpath("container")
cd("container")

pkgs = [
    (name = "ArgParse", version = "1.1.0")
    (name = "BioAlignments", version = "2.0.0")
    (name = "BioSequences", version = "2.0.5")
    (name = "CodecZlib", version = "0.7.0")
    (name = "CSV", version = "0.7.7")
    (name = "DataFrames", version = "0.21.8")
    (name = "DataStructures", version = "0.17.20")
    (name = "Edlib", version = "0.1.1")
    (name = "JSON", version = "0.21.1")
    (name = "FASTX", version = "1.1.3")
    (name = "Query", version = "1.0.0")
    (name = "XAM", version = "0.2.6")
]
no_test = [x.name for x in pkgs]

julia_version = v"1.5.2"
parent_image = "ubuntu:focal"
apt_override = String[
    # absolutely necessary
    "build-essential",
    "ca-certificates",
    "curl",
    "gpg",
    "gpg-agent",
    "locales",
    "wget",
    # additional
    "pigz",
]

# create Dockerfile
SimpleContainerGenerator.create_dockerfile(
    pkgs;
    output_directory = pwd(),
    no_test = no_test,
    julia_version = julia_version,
    parent_image = parent_image,
    override_default_apt = apt_override
)


# build image using Podman, "transform" to Singularity and upload to library
local_image_tag = "dialvarezs/10x_consensus"
library_image_url = "dialvarezs/default/10x-consensus:latest"
basename = split(local_image_tag, "/")[2]

run(`podman build -t $(local_image_tag) .`)
run(`podman save --format oci-archive --output $(basename).tar $(local_image_tag)`)
run(`singularity build $(basename).sif oci-archive://$(basename).tar`)
run(`singularity sign --keyidx 0 $(basename).sif`)
run(`singularity push $(basename).sif library://$(library_image_url)`)
run(`rm -fr $(basename).tar Dockerfile simplecontainergenerator_container_files`)
