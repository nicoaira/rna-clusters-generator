#!/bin/bash

# Define the download URL for rnartistcore
JAR_URL="https://github.com/fjossinet/RNArtistCore/releases/download/0.4.6-SNAPSHOT/rnartistcore-0.4.6-SNAPSHOT-jar-with-dependencies.jar"

# Define the installation directory inside the home directory
INSTALL_DIR="$HOME/.local/bin"
JAR_NAME="rnartistcore.jar"
SCRIPT_NAME="rnartistcore"

# Ensure ~/.local/bin exists
mkdir -p "$INSTALL_DIR"

# Function to check and install Java
install_java() {
    if ! command -v java &> /dev/null; then
        echo "Java is not installed. Installing OpenJDK..."
        sudo apt-get update
        sudo apt-get install -y default-jre
    else
        echo "Java is already installed."
    fi
}

# Function to install ViennaRNA
install_vienna_rna() {
    if ! dpkg -s vienna-rna &> /dev/null; then
        echo "Installing ViennaRNA..."
        sudo apt-get update
        sudo apt-get install -y vienna-rna
    else
        echo "ViennaRNA is already installed."
    fi
}

# Function to install Rust using rustup
install_rust() {
    if ! command -v rustc &> /dev/null; then
        echo "Rust is not installed. Installing Rust using rustup..."
        curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
        export PATH="$HOME/.cargo/bin:$PATH"
    else
        echo "Rust is already installed."
    fi
}

# Install Java, ViennaRNA, and Rust
install_java
install_vienna_rna
install_rust

# Check if wget or curl is available for downloading
if command -v wget &> /dev/null; then
    DOWNLOAD_CMD="wget -q --show-progress -O"
elif command -v curl &> /dev/null; then
    DOWNLOAD_CMD="curl -L -o"
else
    echo "Error: Neither wget nor curl is installed. Please install one and try again."
    exit 1
fi

# Download the JAR file
echo "Downloading rnartistcore JAR..."
$DOWNLOAD_CMD "$INSTALL_DIR/$JAR_NAME" "$JAR_URL"

# Verify if the file was downloaded successfully
if [ ! -f "$INSTALL_DIR/$JAR_NAME" ]; then
    echo "Error: Failed to download the JAR file."
    exit 1
fi

# Create the executable wrapper script
echo "Creating the rnartistcore executable..."
echo '#!/bin/bash' > "$INSTALL_DIR/$SCRIPT_NAME"
echo "java -jar $INSTALL_DIR/$JAR_NAME \"\$@\"" >> "$INSTALL_DIR/$SCRIPT_NAME"

# Make the script executable
chmod +x "$INSTALL_DIR/$SCRIPT_NAME"

# Add ~/.local/bin and Rust to PATH if not already present
if [[ ":$PATH:" != *":$INSTALL_DIR:"* ]]; then
    echo "Adding $INSTALL_DIR to PATH..."
    echo "export PATH=\"$INSTALL_DIR:\$PATH\"" >> "$HOME/.bashrc"
    echo "export PATH=\"$INSTALL_DIR:\$PATH\"" >> "$HOME/.zshrc"  # For zsh users
    export PATH="$INSTALL_DIR:$PATH"
fi

if [[ ":$PATH:" != *"$HOME/.cargo/bin"* ]]; then
    echo "Adding Rust to PATH..."
    echo "export PATH=\"$HOME/.cargo/bin:\$PATH\"" >> "$HOME/.bashrc"
    echo "export PATH=\"$HOME/.cargo/bin:\$PATH\"" >> "$HOME/.zshrc"
    export PATH="$HOME/.cargo/bin:$PATH"
fi

# Verify installation
if command -v rnartistcore &> /dev/null; then
    echo "Installation successful! You can now run rnartistcore from anywhere."
    echo "Example usage:"
    echo "  rnartistcore /path/to/your/script"
    echo "Restart your terminal or run 'source ~/.bashrc' (or 'source ~/.zshrc' for zsh) if the command is not found."
else
    echo "Error: rnartistcore was not installed properly."
    exit 1
fi
