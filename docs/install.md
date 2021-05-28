
# Setup the linux system

The full instructions at [this link](https://docs.microsoft.com/en-us/windows/wsl/install-win10)

Here are the commands in the link above:
1. __Open PowerShell as Administrator.__
2. Enable Windows Subsystem for Linux:

`dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart`

3. Enable Virtual Machine feature:

`dism.exe /online /enable-feature /featurename:VirtualMachinePlatform /all /norestart`

Here's what the last two commands output looks like:

![wsl_commands](/docs/images/wsl_commands.png)

4. Restart.
5. Download Ubuntu from the Microsoft Store. Just click install.

![store_1](/docs/images/store_1.png)

6. Once download click launch.

![store_2](/docs/images/store_2.png)

7. One the first launch you have to create a username and password for the linux machine. Press enter after you type each input. The password does not appear when you type. It will look similar to below.

![first_launch](/docs/images/first_launch.png)

You now have an Ubuntu distribution install.

You can run it by pressing ⊞ Windows key , typing "Ubuntu" , and running the Ubuntu App.

Extra info:

* The to linux filesystem can be found at a path similar to this : `C:\Users\username\AppData\Local\Packages\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\LocalState\rootfs`. You can use the File Explorer to view the contents here. You can copy from here (I have had no problems yet). __Do not add files or folders here with File Explorer.__ Use `cp` or `mv` to copy or move files to the the linux filesystem. Another way is to not move the files and to access the drives directly located at `mnt/<drive_letter>/`. For instance, the C drive is `mnt/c/`.
* You can uninstall by going to __Settings__ > __Apps__, finding Ubuntu and clicking __Uninstall__.

## Pasting in a terminal

Pasting can be done by __right-clicking__ in the terminal. If that does not work you can add copy/paste to the linux terminal by editing the terminal properties.

![edit_wsl_properties](/docs/images/edit_wsl_properties.png)

Select _QuickEdit_ Mode and _Use Ctrl+Shift+C/V as Copy/Paste_. _QuickEdit_ allows for the right-click pasting trick.

![edit_wsl_properties_2](/docs/images/edit_wsl_properties_2.png)

__Ctrl+Shift+C__ and __Ctrl+Shift+V__ will also be copy and paste.

# Install git, python, pip, java, and bowtie

First, update the apt package manager and get current package versions:
```
sudo apt update
```

You can install everything at once with this command:
```
sudo apt install git-all python3 python3-pip openjdk-11-jdk maven bowtie
```
If there are issues, each program can be install separately as shown below.

## Install git

```
sudo apt install git-all
```

Verify with:
```
git --version
```

Here is the help message from `git --version`:

![verify_git](/docs/images/verify_git.png)

## Install python3 and pip

```
sudo apt install python3 python3-pip
```

Verify with:
```
python3 -V
```

Here is the help message from `python3 -V`:

![verify_python](/docs/images/verify_python.png)

## Install java

Trimmomatic uses java. Simply install it with the apt package manager:
```
sudo apt install openjdk-11-jdk maven
```

Verify with:
```
java
```

Here is the help message from `java`:

![verify_java](/docs/images/verify_java.png)


## Install bowtie

Bowtie exists in the apt package manager:
```
sudo apt install bowtie
```

Verify with:
```
bowtie
bowtie-build
bowtie2-inspect
```

Here is the help message from `bowtie2-inspect`:

![verify_bowtie](/docs/images/verify_bowtie.png)
