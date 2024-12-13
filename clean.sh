#!/bin/sh
current_directory="$(cd "$(dirname "$0")" && pwd)"
project_root_dir="$current_directory"
rm -rf "${project_root_dir}/third_party"
rm -rf "${project_root_dir}/build"
rm -rf "${project_root_dir}/install"
rm -rf "${project_root_dir}/results"
rm -rf "${project_root_dir}/docs"

echo "declare success -- super duper clean now :P"