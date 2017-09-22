# Matrix
by Paul Herz

I developed this library off of Dr. Daniel Reynold's implementation. It has been approved as a drop-in replacement for that library in the course of the High Performance Scientific Computing class.

## Use
To use the included `Matrix` and `Vector` classes in your project, simply drop the contents of the `src` folder here into your project. You may want to separate this code into a `lib` folder to denote that it is not your project code if using it in HPSC.

### Including into a traditional project (GCC, Makefiles)
1. Copy the contents of `src` in this repo, including all headers, into a location in your project.
2. Add `Matrix.cpp` and `Vector.cpp` to your build list (e.g. `g++ yourfile.cpp otherfile.cpp Matrix.cpp Vector.cpp`, including the relative path if they are in a different folder than your code.
3. If you've placed this code in a different folder than where you are running the compiler, add the include header folder flag (`-I /lib` if the files are in `lib`).
4. You will need to use at least the C++11 standard, specifically the GNU extension, using the flag `-std=gnu++11`. Your project may not work if it is using an older standard or a non-GNU extension standard. If you have trouble with a *newer* standard, file an issue and I will respond promptly.

### Including into an Xcode Project (Xcode 9.0)
This example starts from absolute scratch in Xcode 9.0. This version changed the process of including files to be a little more complicated and restrictive. As such, it is difficult to actually store this code in a separate folder from your own. However, we can still use Xcode's Groups (pseudo-folders in the Xcode sidebar) to *visually* group them apart from your own code.

1. Open Xcode to the welcome screen, and click **Create a New Xcode Project**. If you don't get a welcome screen when you open Xcode, go to `File > New > Project...` in the toolbar, or press `Shift + Cmd + N`.
2. Select **macOS** from the tabs at the top of the New Project wizard, then select **Command Line Tool** under the **Application** category (at the top of the list). Click **Next**.
3. Enter a product name, preferably without spaces. Select **C++** as your language. Click **Next**.
4. Pick a location to save your project to. Be wary of the **Create Git repository on my Mac**. You may or may not want this checked. Xcode will create a folder named after your project here, so there's no need to do so manually.
5. Open up a Finder window next to Xcode. In that window, browse to where you've downloaded/cloned this code, specifically to the `src` folder. Select all of the *contents* of this folder (not the folder itself) and *drag* the files into Xcode, dropping them directly above or below the `main.cpp` file. 
6. In the popup, select **Copy items if needed** and **Create groups**. Click **Finish**.
7. These files may clutter your future or existing code. You can split them off into a group on the sidebar, but it will not change their actual location in the folder. Select all the files belonging to the Matrix library by clicking the first file, holding `Shift`, and clicking the last. Right click, and select **New Group from Selection**. Name it something like `lib` or `matrix` to identify it as separate. Now all of the Matrix library files will be bundled into a folder and separated. If you don't want this group *inside* of your `ProjectName` group, you can carefully drag it above the `ProjectName` group on the sidebar in Xcode and the group as well as the folder will be moved to a *sibling* of the `ProjectName` group.
8. If you look in Finder and it appears that Xcode has *copied*, rather than *moved* the Matrix files in your project (there is one copy of all the files inside /ProjectName/ProjectName and another in /ProjectName/ProjectName/GroupName), you will need to manually delete the old copies (in /ProjectName/ProjectName). This is a behavioral issue in Xcode 9.0.
9. Select your project from the file sidebar in Xcode (it's the blueprint icon, at the top usually). Select the **Build Settings** tab. Search `C++ Language Dialect` in the search bar below the tab bar.
10. In the resulting setting that appears, make sure the selection is *at least* **GNU++11 [-std=gnu++11]**. This library depends on at least C++11, and has been tested with only the C++11 standard with GNU extensions. Further testing is to be done.
11. Select the **Build Phases** tab now under the same project, and open the **Compile Sources** dropdown by clicking on the arrow next to it. Make sure `Matrix.cpp` and `Vector.cpp` are added. If they are not, click the plus and add them.
12. You should be able to add the following to your `main.cpp` file or any other file in the project:
```
#include "Matrix.h"
```
Xcode should be able to find the header files.
