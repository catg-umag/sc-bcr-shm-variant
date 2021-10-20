#!/usr/bin/env julia
using SimpleContainerGenerator

mkpath("container")
cd("container")

pkgs = [
    (name = "ArgMacros", version = "0.2.4")
    (name = "ArgParse", version = "1.1.4")
    (name = "BioAlignments", version = "2.0.0")
    (name = "BioSequences", version = "2.0.5")
    (name = "CodecZlib", version = "0.7.0")
    (name = "CSV", version = "0.8.5")
    (name = "DataStructures", version = "0.18.10")
    (name = "Edlib", version = "0.1.1")
    (name = "FASTX", version = "1.2.0")
    (name = "JSON3", version = "1.9.2")
    (name = "Query", version = "1.0.0")
    (name = "StructTypes", version = "1.8.1")
    (name = "XAM", version = "0.2.7")
]
no_test = [x.name for x in pkgs]

julia_version = v"1.6.3"
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
