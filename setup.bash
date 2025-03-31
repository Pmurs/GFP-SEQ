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

# Make sure Python 3.9 is installed
sudo apt install -y python3.9 python3.9-venv

# Create a new venv using Python 3.9
sudo make