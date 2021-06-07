#!/usr/bin/env julia
using SimpleContainerGenerator

mkpath("container")
cd("container")

pkgs = [
    (name = "ArgMacros", version = "0.2.3")
    (name = "ArgParse", version = "1.1.4")
    (name = "BioAlignments", version = "2.0.0")
    (name = "BioSequences", version = "2.0.5")
    (name = "CodecZlib", version = "0.7.0")
    (name = "CSV", version = "0.8.5")
    (name = "DataStructures", version = "0.18.9")
    (name = "Edlib", version = "0.1.1")
    (name = "FASTX", version = "1.1.3")
    (name = "JSON3", version = "1.8.2")
    (name = "Query", version = "1.0.0")
    (name = "StructTypes", version = "1.7.2")
    (name = "XAM", version = "0.2.7")
]
no_test = [x.name for x in pkgs]

julia_version = v"1.6.1"
parent_image = "ubuntu:20.04"
apt_override = String[
    # absolutely necessary
    "build-essential",
    "ca-certificates",
    "curl",
    "gpg",
    "gpg-agent",
    "locales",
    "procps",
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
    override_default_apt = apt_override,
)


# build image using Podman, "transform" to Singularity and upload to library
local_image_tag = "dialvarezs/10x_consensus"
library_image_url = "dialvarezs/default/10x-consensus:latest"
basename = split(local_image_tag, "/")[2]

run(`podman build -t $(local_image_tag) .`)
run(`podman save --format oci-archive --output $(basename).tar $(local_image_tag)`)
run(`singularity build -F $(basename).sif oci-archive://$(basename).tar`)
run(`rm -fr $(basename).tar Dockerfile simplecontainergenerator_container_files`)
run(`singularity sign $(basename).sif`)
run(`singularity push $(basename).sif library://$(library_image_url)`)
