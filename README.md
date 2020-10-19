
This is a guide for installing and running bowtie2 and localmapper. I made this for Windows 10 users who don't want to become computer experts, but still need use bowtie2 efficiently. The main goal is to make sure you are running the right software with the appropriate inputs.

bowtie2: [https://github.com/BenLangmead/bowtie2](https://github.com/BenLangmead/bowtie2)
---
localmapper: [https://github.com/DNAReplicationLab/localMapper](https://github.com/DNAReplicationLab/localMapper)
---

# Contents

- [Setup the linux system](#setup-the-linux-system)
- [Install bowtie2](#install-bowtie2)
- [Install localmapper](#install-localmapper)
- [Example workflows](#example-workflows)
- [Extra Linux commands](#extra-linux-commands)

# Setup the linux system

The full instructions at [this link](https://docs.microsoft.com/en-us/windows/wsl/install-win10)

Here are the commands in the link above:
1. __Open PowerShell as Administrator.__
2. Enable Windows Subsystem for Linux:

`dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart`

3. Enable Virtual Machine feature:

`dism.exe /online /enable-feature /featurename:VirtualMachinePlatform /all /norestart`

Here's what the last two commands output looks like:

![wsl_commands](/images/wsl_commands.png)

4. Restart.
5. Download Ubuntu from the Microsoft Store. Just click install.

![store_1](/images/store_1.png)

6. Once download click launch.

![store_2](/images/store_2.png)

7. One the first launch you have to create a username and password for the linux machine. Press enter after you type each input. The password does not appear when you type. It will look similar to below.

![first_launch](/images/first_launch.png)

You now have an Ubuntu distribution install.

You can run it by pressing ⊞ Windows key , typing "Ubuntu" , and running the Ubuntu App.

Extra info:

* The to linux filesystem can be found at a path similar to this : `C:\Users\username\AppData\Local\Packages\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\LocalState\rootfs`. You can use the File Explorer to view the contents here. You can copy from here (I have had no problems yet). __Do not add files or folders here with File Explorer.__ Use `cp` or `mv` to copy or move files to the the linux filesystem. Another way is to not move the files and to access the drives directly located at `mnt/drive_letter/`. For instance, the C drive is `mnt/c/`.
* You can uninstall by going to __Settings__ > __Apps__, finding Ubuntu and clicking __Uninstall__.

## Pasting in a terminal

Pasting can be done by __right-clicking__ in the terminal. If that does not work you can add copy/paste to the linux terminal by editing the terminal properties.

![edit_wsl_properties](/images/edit_wsl_properties.png)

Select _QuickEdit_ Mode and _Use Ctrl+Shift+C/V as Copy/Paste_. _QuickEdit_ allows for the right-click pasting trick.

![edit_wsl_properties_2](/images/edit_wsl_properties_2.png)

__Ctrl+Shift+C__ and __Ctrl+Shift+V__ will also be copy and paste.

## Setup the project folders on Windows

I store my whole project in one folder on my C drive. I put at at this path :`C:\projects\bowtie2`

Then, I make the following folders. They are all empty now.

![project_folder](/images/project_folder.png)

# Install bowtie2

There are a few ways to get bowtie2. I use the first method.

__Method 1 - apt package manager__

1. Start the Ubuntu app. It will open a terminal.

![open_terminal](/images/open_terminal.png)

2. Run `sudo apt install bowtie2`. You will need to enter your password. The terminal spaces will remain empty as you type your password.

![apt_bowtie2](/images/apt_bowtie2.png)

3. Done

__Method 2 - get compiled binaries__

No instructions yet.

## Verify bowtie2 install

The binaries are stored in /bin or /usr/bin and they called:
```
bowtie2
bowtie2-build
bowtie2-inspect
```
Run all three of the commands listed above in any directory of the linux terminal. All commands should display help messages.

Here is the help message from `bowtie2-inspect`:

![message_1](/images/message_1.png)


# Install localmapper

We need to get the localmapper scripts from the github repo into our project folder. There are 2 methods.

__Method 1 - You have git installed on windows__

A short cut for opening a terminal at a folder is to __Shift+Right Click__ in the explorer window and click __Open Powershell window here__.

![open_cmd_from_folder](/images/open_cmd_from_folder.png)

Run `git clone https://github.com/DNAReplicationLab/localMapper.git` in your project directory.

![git_clone](/images/git_clone.png)

__Method 2 - You do not have git installed__

1. Download the repo as zip.

The source repo can be found [here](https://github.com/DNAReplicationLab/localMapper)

![download_repo](/images/download_repo.png)

2. Extract and copy it to your project folder (mine is `C:\projects\bowtie2`).

Both methods result in a folder called `localMapper` in your porject folder and it will have identical contents to the repo.

Get additional packages for localmapper: `sudo apt-get install bowtie2 samtools bedtools picard-tools`

**Not finshed

## Verify localmapper install

Note : always run the localmapper from the directory where you want the outputs to go. I always run mine in `/home/user/bowtie2/`


# Example workflows

## Example 1 - Running the examples in the bowtie2 guide
1. Downloading the example files



## Example 2 - Getting started with localmapper


## Example 3 - Lab data


# Extra Linux commands

## Stopping programs

There are a few ways to stop programs running in the linux terminal. Try them in this order;
1. click in the terminal to activate the window, then spam _Ctrl+C_ a bunch to see what happens. Most programs will exit and display an inturrept message.
2. If _Crtl+C_ does not work, then you can just exit out of the terminal window.
3. If the window is frozen, open the tsak manager (_Ctrl+Shift+Esc_), __right-click__ on Ubuntu, click __End task__.

## Printing the path

The path variable is a bunch of folder locations where executable programs can be found. You can view it with `echo $PATH`, but I think that is hard to read. Each location is separated by a colon.

I found this on stackoverflow (I don't have the link) : `echo $PATH | tr ":" "\n" | nl`

