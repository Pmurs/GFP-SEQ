#!/bin/bash

# Fast install of nanopore alignment tools using apt
# Assumes you're in a Python 3 virtual environment

echo "Updating apt package list..."
sudo apt update

echo "Installing core tools via apt..."
sudo apt install -y \
  filtlong \
  minimap2 \
  samtools \
  fastqc
medaka

echo "✅ All tools installed!"
